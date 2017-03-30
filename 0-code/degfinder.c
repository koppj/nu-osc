/* ---------------------------------------------------------------------------- */
/* Finding degeneracies in neutrino oscillations                                */
/* ---------------------------------------------------------------------------- */
/* Author: Joachim Kopp (partly based on ideas by P. Huber and W. Winter)       */
/* ---------------------------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <globes/globes.h>
#include "nu.h"
#include "snu.h"
#ifdef NU_USE_MONTECUBES
  #include <montecubes/montecubes.h>
#endif


// ----------------------------------------------------------------------------
double ChiNPWrapper(glb_params base_values, int hierarchy, glb_params fit_values)
// ----------------------------------------------------------------------------
// Perform minimization for a certain set of test values (base_values) and
// return the resulting fit values (fit_values). The user has to specify the
// desired hierarchy (HIERARCHY_NORMAL or HIERARCHY_INVERTED).
// ----------------------------------------------------------------------------
{
  double result;
  glb_params tv     = glbAllocParams();      /* Test values             */
  glb_params cv     = glbAllocParams();      /* Central values          */
  glb_params old_cv = glbAllocParams();      /* Previous central values */
  glbCopyParams(base_values, tv);
  glbGetCentralValues(cv);
  glbGetCentralValues(old_cv);

  // Note: Converting the dm31^2 for the normal hierarchy into a dm31^2 for the
  // inverted hierarchy that would give similar oscillation probabilities has bee
  // discussed in hep-ph/0503079 (see also hep-ph/0509359 and 0812.1879).
  // While a suitable transformation can be found approximately, it is different
  // for different oscillation channels. In all cases, it is similar to 2*s12^2*dm21^2
  // or 2*c12^2*dm21^2, which, very roughly, is of O(dm21^2), so that's what
  // we use here.
  double dm31_tv = glbGetOscParams(tv,GLB_DM_31);
  double dm31_cv = glbGetOscParams(cv,GLB_DM_31);
  if (hierarchy == HIERARCHY_NORMAL)
  {
    if (dm31_tv < 0)  glbSetOscParams(tv, -dm31_tv + glbGetOscParams(tv,GLB_DM_21), GLB_DM_31);
    if (dm31_cv < 0)  glbSetOscParams(cv, -dm31_cv + glbGetOscParams(cv,GLB_DM_21), GLB_DM_31);
    glbSetCentralValues(cv);

    // Before calling glbChiNP, check if GLoBES will accept the chosen
    // oscillation parameters (it may not, for instance if the chosen values of
    // s22thmue and Um4 are inconsistent)
    if (glbSetOscillationParameters(tv) == 0)
      result = glbChiNP(tv, fit_values, GLB_ALL);
    else
      result = 1.e21;
  }
  else if (hierarchy == HIERARCHY_INVERTED)
  {
    if (dm31_tv > 0)  glbSetOscParams(tv, -dm31_tv + glbGetOscParams(tv,GLB_DM_21), GLB_DM_31);
    if (dm31_cv > 0)  glbSetOscParams(cv, -dm31_cv + glbGetOscParams(cv,GLB_DM_21), GLB_DM_31);
    glbSetCentralValues(cv);
    if (glbSetOscillationParameters(tv) == 0)
      result = glbChiNP(tv, fit_values, GLB_ALL);
    else
      result = 1.e22;
  }
  else
  {
    fprintf(stderr, "ChiNPWrapper: Please specify HIERARCHY_NORMAL or HIERARCHY_INVERTED!\n");
    result = -1.0;
  }

  glbSetCentralValues(old_cv);
  glbFreeParams(old_cv);
  glbFreeParams(cv);
  glbFreeParams(tv);

  return result;
}


/* ---------------------------------------------------------------------------- */
int degfinder(const glb_params base_values, const int n_prescan_params,
      const int *prescan_params, const double *prescan_min,
      const double *prescan_max, const int *prescan_steps,
      const glb_projection prescan_proj, const glb_projection fit_proj,
      int *n_deg, glb_params *deg_pos, double *deg_chi2, const unsigned long flags,
      const unsigned long *prescan_flags, const char *output_file)
/* ---------------------------------------------------------------------------- */
/* Input parameters:                                                            */
/*   base_values: The oscillation parameters                                    */
/*   n_prescan_params: Number of parameters to perform prescan on               */
/*   prescan_params: Indices of the parameters for the prescan                  */
/*   prescan_min, prescan_max, prescan_steps: Minimum/Maximum values and        */
/*     numbers of steps for the prescan                                         */
/*   prescan_proj, fit_proj: Projections for the prescan and for the final fit  */
/*   flags: A combination of the DEG_XXX flags                                  */
/*   prescan_flags: Options for each of the parameters that are scanned         */
/*     (e.g. DEG_LOGSCALE, DEG_PLUS_MINUS, DEG_S22)                             */
/*   output_file: path and base name for output files (in MCMC mode)            */
/*     if NULL, file name will be chosen automatically                          */
/* Output parameters:                                                           */
/*   n_deg: Input:  Max. number of degenerate solutions to return               */
/*          Output: Number of degenerate solutions found                        */
/*   deg_pos: Positions of degeneracies in parameter space                      */
/*   deg_chi2: chi^2 values of the degenerate solutions                         */
/* ---------------------------------------------------------------------------- */
{
  /* Copy input parameters to private data structures */
  long private_flags = flags;
  int n_p_params = n_prescan_params;
  int p_params[n_p_params];
  double p_max[n_p_params], p_min[n_p_params];
  int p_steps[n_p_params];
  int p_flags[n_p_params];
  for (int i=0; i < n_p_params; i++)
  {
    p_params[i] = prescan_params[i];
    p_min[i]    = prescan_min[i];
    p_max[i]    = prescan_max[i];
    p_steps[i]  = prescan_steps[i];
    if (!prescan_flags)
      p_flags[i] = 0;
    else
      p_flags[i] = prescan_flags[i];
  }
  glb_projection private_prescan_proj = glbAllocProjection();
  glb_projection private_fit_proj     = glbAllocProjection();
  glbCopyProjection(prescan_proj, private_prescan_proj);
  glbCopyProjection(fit_proj, private_fit_proj);
  
  /* If no explicit flags given, determine automatically which parameters
   * should be scanned on a log scale */
  if (!prescan_flags)
  {
    for (int i=0; i < n_p_params; i++)
    {
      if (p_params[i] == GLB_THETA_13)
        p_flags[i] |= (DEG_LOGSCALE | DEG_S22);
      else if (strstr(snu_param_strings[p_params[i]], "ARG")   == NULL ||
               strstr(snu_param_strings[p_params[i]], "DELTA") == NULL)
      {
        p_flags[i] |= (DEG_LOGSCALE | DEG_PLUS_MINUS);

        /* If phase of this parameter is explicitly included -> OK. Otherwise, use
         * DEG_PLUS_MINUS flag to at least include the sign ambiguity */
        for (int j=0; j < n_p_params; j++)
        {
          if (p_params[j] == p_params[i]+1  &&  strstr(snu_param_strings[p_params[j]], "ARG") != NULL)
            p_flags[i] &= (~DEG_PLUS_MINUS);
        }
      }
    } /* for(i) */
  }

  /* For runs without systematics, switch them off now */
  int old_sys_state[glb_num_of_exps][128];
  for (int j=0; j < glb_num_of_exps; j++)
    for (int k=0; k < glbGetNumberOfRules(j); k++)
      old_sys_state[j][k] = glbGetSysOnOffState(j, k);
  if (private_flags & DEG_NO_SYS)
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);

  /* For runs without degeneracies, omit prescan and wrong hierarchy solutions */
  if (private_flags & DEG_NO_DEG)
  {
    n_p_params = 0;
    if (glbGetOscParams(base_values, GLB_DM_31) > 0)
      private_flags |= DEG_NO_IH;
    else
      private_flags |= DEG_NO_NH;
  }
  
  /* For runs without non-standard degeneracies, remove all NSI parameters from
   * parameter list */
  if (private_flags & DEG_ONLY_STD_DEG)
  {
    for (int i=0; i < n_p_params; i++)
    {
      if (p_params[i] >= 6)
      {
        for (int j=i; j < n_p_params-1; j++)
        {
          p_params[j] = p_params[j+1];
          p_min[j]    = p_min[j+1];
          p_max[j]    = p_max[j+1];
        }
        n_p_params--;
      }
    }
  }

  /* For runs without correlations, switch them off */
  if (private_flags & DEG_NO_CORR)
  {
    for (int i=0; i < glbGetNumOfOscParams(); i++)
    {
      glbSetProjectionFlag(private_prescan_proj, GLB_FIXED, i);
      glbSetProjectionFlag(private_fit_proj, GLB_FIXED, i);
    }
    glbSetDensityProjectionFlag(private_prescan_proj, GLB_FIXED, GLB_ALL);
    glbSetDensityProjectionFlag(private_fit_proj, GLB_FIXED, GLB_ALL);
  }


  /* Preparations for prescan */
  /* ------------------------ */
  
  /* For cyclic parameters (complex phases), make sure the upper boundary
   * of the scan region is omitted if it would imply double counting */
  for (int i=0; i < n_p_params; i++)
  {
    p_steps[i] = prescan_steps[i];
    p_min[i]   = prescan_min[i];
    if ( p_params[i] == GLB_DELTA_CP ||
         strstr(snu_param_strings[p_params[i]], "ARG") != NULL ||
         strstr(snu_param_strings[p_params[i]], "DELTA") != NULL )
    {
      if (fabs(prescan_min[i]) < 1e-10 && fabs(prescan_max[i] - 2*M_PI) < 1e-10)
      {
        p_max[i]    = prescan_max[i] - (prescan_max[i] - prescan_min[i])/prescan_steps[i];
        p_steps[i] -= 1;
      }
      else
        p_max[i] = prescan_max[i];
    }
    else
      p_max[i] = prescan_max[i];
  }

  /* If we are asked to consider positive and negative parameter
   * values, we have to use twice as many prescan points */
  for (int i=0; i < n_p_params; i++)
    if (p_flags[i] & DEG_PLUS_MINUS)
      p_steps[i] = 2*p_steps[i] + 1;
 
  /* Compute number of points in prescan grid */
  unsigned long n_prescan_points = 1;
  for (int i=0; i < n_p_params; i++)
    n_prescan_points *= p_steps[i] + 1;

  /* Function for converting one-dimensional indices to a multi-dimensional
   * index for (d+1)-th dimension */
  int convert_index(int j, int d)
  {
    int k = n_prescan_points;
    for (int i=0; i <= d; i++)
    {
      j %= k;
      k /= p_steps[i] + 1;
    }
    j /= k;
    return j;
  }

  if (debug_level > 1)
    printf("#   Using %lu prescan points.\n", n_prescan_points);
  
  /* Create data structures */
  double chi2_table_NH[n_prescan_points];
  double chi2_table_IH[n_prescan_points];
  glb_params Fit_NH = glbAllocParams();
  glb_params Fit_IH = glbAllocParams();
  glb_params mcb_conv_crit = glbAllocParams();
  glb_params mcb_steps = glbAllocParams();
  glb_params param_table_NH[n_prescan_points];
  glb_params param_table_IH[n_prescan_points];
  for (int j=0; j < n_prescan_points; j++)
  {
    param_table_NH[j] = glbAllocParams();
    param_table_IH[j] = glbAllocParams();    
  }
  memset(chi2_table_NH, 0.0, sizeof(chi2_table_NH[0]) * n_prescan_points);
  memset(chi2_table_IH, 0.0, sizeof(chi2_table_NH[0]) * n_prescan_points);

  if (debug_level > 1)
  {
    printf("#   Parameters used in prescan: ");
    for (int i=0; i < n_p_params; i++)
    {
      printf("%s ", snu_param_strings[p_params[i]]);
      if (p_flags[i] & DEG_LOGSCALE)
        printf("(log) ");
      if (p_flags[i] & DEG_PLUS_MINUS)
        printf("(+-) ");
      if (p_flags[i] & DEG_S22)
        printf("(s22) ");
    }
    printf("\n");
  }

  
  /* Prescan */
  /* ------- */
 
  glbSetProjection(private_prescan_proj); 
  glbCopyParams(base_values, Fit_NH);
  glbCopyParams(base_values, Fit_IH);
  for (int j=0; j < n_prescan_points; j++)
  {
    /* Compute parameter values at current grid point */
    double prescan_test_values[n_p_params];
    for (int i=0; i < n_p_params; i++)
    {
      double x;
      if (p_flags[i] & DEG_LOGSCALE)        /* Use log scale for this param? */
      {
        if (p_flags[i] & DEG_PLUS_MINUS)
        {
          int k = convert_index(j,i);
          if (k <= (p_steps[i]-1) / 2)
            x =  -POW10( p_max[i] - k * (p_max[i]-p_min[i])/(p_steps[i]/2) );
          else
            x =   POW10( p_min[i] + (k - (p_steps[i]+1)/2) * (p_max[i]-p_min[i])/(p_steps[i]/2) );
        }
        else
          x = POW10( p_min[i] + convert_index(j,i) * (p_max[i]-p_min[i])/p_steps[i] );
      }
      else                                  /* Otherwise: linear distribution */
      {
        if (p_flags[i] & DEG_PLUS_MINUS)
        {
          int k = convert_index(j,i);
          if (k <= (p_steps[i]-1) / 2)
            x = -( p_max[i] - k * (p_max[i]-p_min[i])/(p_steps[i]/2) );
          else
            x =    p_min[i] + (k - (p_steps[i]+1)/2) * (p_max[i]-p_min[i])/(p_steps[i]/2);
        }
        else
          x = p_min[i] + convert_index(j,i) * (p_max[i]-p_min[i])/p_steps[i];

        prescan_test_values[i] = x;
      }
      if (p_flags[i] & DEG_S22)
        prescan_test_values[i] = SGN(x) * asin(sqrt(fabs(x)))/2.0;
      else
        prescan_test_values[i] = x;

      glbSetOscParams(Fit_NH, prescan_test_values[i], p_params[i]);
      glbSetOscParams(Fit_IH, prescan_test_values[i], p_params[i]);
    }

    /* Compute chi^2 _without_ systematics for both hierarchies */
    int old_sys_state2[glb_num_of_exps][128];
    for (int j=0; j < glb_num_of_exps; j++)
      for (int k=0; k < glbGetNumberOfRules(j); k++)
        old_sys_state2[j][k] = glbGetSysOnOffState(j, k);
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
    if ( !(private_flags & DEG_NO_NH) )
    {
      chi2_table_NH[j] = ChiNPWrapper(Fit_NH, HIERARCHY_NORMAL, Fit_NH);
      glbCopyParams(Fit_NH, param_table_NH[j]);
    }
    if ( !(private_flags & DEG_NO_IH) )
    {
      chi2_table_IH[j] = ChiNPWrapper(Fit_IH, HIERARCHY_INVERTED, Fit_IH);
      glbCopyParams(Fit_IH, param_table_IH[j]);
    }
    for (int j=0; j < glb_num_of_exps; j++)
      for (int k=0; k < glbGetNumberOfRules(j); k++)
        glbSwitchSystematics(j, k, old_sys_state2[j][k]);

    if (debug_level > 2)
    {
      printf("#   Test params:");
      for (int i=0; i < n_p_params; i++)
        printf(" %8.5g", prescan_test_values[i]);
      printf(", chi2_NH = %12.7g, chi2_IH = %12.7g\n", chi2_table_NH[j], chi2_table_IH[j]);
    }
  }

  
  /* Determine locations of the degeneracies from chi^2 tables and perform fit */
  /* ------------------------------------------------------------------------- */

  glbSetProjection(private_fit_proj);
  int max_deg = *n_deg;
  *n_deg = 0;
  for (int j=0; j < n_prescan_points; j++)
  {
    /* Compute differences to the chi^2 values from neighbouring points to
     * find minima */
    int k = n_prescan_points;
    int min_NH = 1, min_IH = 1;
    for (int i=0; i < n_p_params; i++)
    {
      int k2 = k;
      k /= (p_steps[i] + 1);
      int s = (j+k) / k2 - j / k2;
      int t = (j+k2) / k2 - (j+k2-k) / k2;
                          /* This is 1 if wraparound happens, and 0 otherwise */
      
      /* Treat phases as cyclic prameters */
      if ( p_params[i] == GLB_DELTA_CP ||
           strstr(snu_param_strings[p_params[i]], "ARG") != NULL ||
           strstr(snu_param_strings[p_params[i]], "DELTA") != NULL )
      {
        /* Ignore chi^2 > 1e10 since this usually indicates some error */
        if (chi2_table_NH[j] > 1.e10 ||
            chi2_table_NH[j] - chi2_table_NH[j + k - s*k2] > 0 ||
            chi2_table_NH[j] - chi2_table_NH[j - k + t*k2] > 0)
          min_NH = 0;
        if (chi2_table_IH[j] > 1.e10 ||
            chi2_table_IH[j] - chi2_table_IH[j + k - s*k2] > 0 ||
            chi2_table_IH[j] - chi2_table_IH[j - k + t*k2] > 0)
          min_IH = 0;
      }
      else
      {
        if (chi2_table_NH[j] > 1.e10 ||
            chi2_table_NH[j] - chi2_table_NH[MIN(n_prescan_points - (k-j%k), j+k)] > 0 ||
            chi2_table_NH[j] - chi2_table_NH[MAX(0, j-k)] > 0)
          min_NH = 0;
        if (chi2_table_IH[j] > 1.e10 ||
            chi2_table_IH[j] - chi2_table_IH[MIN(n_prescan_points - (k-j%k), j+k)] > 0 ||
            chi2_table_IH[j] - chi2_table_IH[MAX(0, j-k)] > 0)
          min_IH = 0;
      }
    }
    
    if ( !(private_flags & DEG_NO_NH) && min_NH )  /* We have found a minimum for the normal hierarchy */
    {
      if (debug_level > 1)
      {
        printf("#   Degeneracy NH at: ");
        for (int k=0; k < n_p_params; k++)
          printf("%d (%s = %g)    ", convert_index(j, k), snu_param_strings[p_params[k]],
                 glbGetOscParams(param_table_NH[j], p_params[k]));
        printf("\n");
        printf("#     Prescan: ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(param_table_NH[j], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(param_table_NH[j], p_params[k]));
        printf(", chi2 = %10.5g\n", chi2_table_NH[j]);
      }

      if (*n_deg >= max_deg)
      {
        fprintf(stderr, "degfinder: Too many degeneracies fond (max. is %d).\n", max_deg);
        break;
      }
      if (flags & DEG_MCMC)  // MCMC mode: remember degenerate solution
        glbCopyParams(param_table_NH[j], deg_pos[*n_deg]);
      else
        deg_chi2[*n_deg] = ChiNPWrapper(param_table_NH[j], HIERARCHY_NORMAL, deg_pos[*n_deg]);

      if (debug_level > 1)
      {
        printf("#     Fit:     ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], p_params[k]));
        printf(", chi2 = %10.5g\n", deg_chi2[*n_deg]);
      }
      (*n_deg)++;
    }
    
    if ( !(private_flags & DEG_NO_IH) && min_IH )  /* We have found a minimum for the inverted hierarchy */
    {
      if (debug_level > 1)
      {
        printf("#   Degeneracy IH at: ");
        for (int k=0; k < n_p_params; k++)
          printf("%d (%s = %g)    ", convert_index(j, k), snu_param_strings[p_params[k]],
                 glbGetOscParams(param_table_IH[j], p_params[k]));
        printf("\n");
        printf("#     Prescan: ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(param_table_IH[j], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(param_table_IH[j], p_params[k]));
        printf(", chi2 = %10.5g\n", chi2_table_IH[j]);
      }

      if (*n_deg >= max_deg)
      {
        fprintf(stderr, "degfinder: Too many degeneracies fond (max. is %d).\n", max_deg);
        break;
      }
      if (flags & DEG_MCMC)
        glbCopyParams(param_table_IH[j], deg_pos[*n_deg]);
      else
        deg_chi2[*n_deg] = ChiNPWrapper(param_table_IH[j], HIERARCHY_INVERTED, deg_pos[*n_deg]);

      if (debug_level > 1)
      {
        printf("#     Fit:     ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], p_params[k]));
        printf(", chi2 = %10.5g\n", deg_chi2[*n_deg]);
      }
      (*n_deg)++;
    }
  } // for(j)


  // MCMC mode: tell MonteCUBES about degeneracies, run MCMC
  // -------------------------------------------------------
  if (flags & DEG_MCMC)
  {
#ifdef NU_USE_MONTECUBES
    mcb_setChainNo(4);                             // Number of MCMC chains
    mcb_setBurnNo(MCB_DYNAMIC_BURN);
    mcb_setLengthMax(1e7);                         // Max. length of each chain
//FIXME    mcb_setLengthMin(5000);                        // Min. length of each chain
    mcb_setLengthMin(50);                        // Min. length of each chain
    mcb_setVerbosity(5);
    mcb_addStartPosition(base_values);
    int deg_pos_ind[*n_deg];
    for (int i=0; i < *n_deg; i++)
      deg_pos_ind[i] = i;
    mcb_setDegeneracySteps(deg_pos, deg_pos_ind, *n_deg); // Tell MonteCUBES about degeneracies
//    mcb_setVerbosity(100);//FIXME
//    mcb_seedMersenne(777);//FIXME
    for (int i=0; i < glbGetNumOfOscParams(); i++) // MonteCUBES convergence criteria
//      glbSetOscParams(mcb_conv_crit, 0.025, i);//FIXME
      glbSetOscParams(mcb_conv_crit, 2.0, i);//FIXME
    glbSetDensityParams(mcb_conv_crit, 1.0, GLB_ALL);
    mcb_setConvergenceCriteria(mcb_conv_crit);

    for (int i=0; i < glbGetNumOfOscParams(); i++) // MonteCUBES step size
    {
      char *p = glbGetParamName(i);
      if (strstr(p, "TH") != NULL)
        glbSetOscParams(mcb_steps, M_PI/30., i);
      else if (strstr(p, "DM") != NULL)
        if (glbGetOscParams(base_values, i) < 1e-20)
          glbSetOscParams(mcb_steps, 0.3, i);
        else
          glbSetOscParams(mcb_steps, 0.01*glbGetOscParams(base_values, i), i);
      else if (strstr(p, "ARG") != NULL  ||  strstr(p, "DELTA") != NULL)
        glbSetOscParams(mcb_steps, M_PI/5., i);
      else
        glbSetOscParams(mcb_steps, 0.001, i);
      free(p);
    }
    glbSetDensityParams(mcb_steps, 0.001, GLB_ALL);
    mcb_setStepSizes(mcb_steps);

    char s[FILENAME_MAX];
    if (!output_file)
      snprintf(s, FILENAME_MAX, "nu-mcmc");
    else
    {
      if (strlen(output_file) > FILENAME_MAX-1)
        fprintf(stderr, "degfinder: name of output file too long.\n");
      strncpy(s, output_file, FILENAME_MAX);
      s[FILENAME_MAX-1] = '\0';
    }

    mcb_full_MCMC(s, GLB_ALL, GLB_ALL, MCB_FULL_MARGIN);
//    mcb_MCMC("/temp/jkopp/sterile-nu/mcmc/mcb-test-oldmcmc.dat", GLB_ALL, GLB_ALL);
#else
    fprintf(stderr, "degfinder: MonteCUBES support not compiled.\n");
#endif
  }
  
  
  /* Clean up */
  if (private_flags & DEG_NO_SYS)
  {
    for (int j=0; j < glb_num_of_exps; j++)
      for (int k=0; k < glbGetNumberOfRules(j); k++)
        glbSwitchSystematics(j, k, old_sys_state[j][k]);
  }
  for (int j=0; j < n_prescan_points; j++)
  {
    if (param_table_NH[j] != NULL)
      glbFreeParams(param_table_NH[j]);
   if (param_table_IH[j] != NULL)
      glbFreeParams(param_table_IH[j]);
  }
  glbFreeParams(mcb_steps);
  glbFreeParams(mcb_conv_crit);
  glbFreeParams(Fit_NH);
  glbFreeParams(Fit_IH);
  glbFreeProjection(private_prescan_proj);
  glbFreeProjection(private_fit_proj);
  
  return 0;
}


