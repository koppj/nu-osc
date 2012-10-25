

// FIXME NOT WORKING! 



/***************************************************************************
 * Non-standard interactions in a neutrino factory with silver channel     *
 ***************************************************************************
 * Author: Joachim Kopp, Toshihiko Ota, Walter Winter                      *
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <globes/globes.h>   /* GLoBES library */
#define HIERARCHY_NORMAL     1
#define HIERARCHY_INVERTED  -1
#define MAX_DEG    100    /* Maximum number of degeneracies to expect */
#define MAX_PARAMS  20

/*************************************************************************** 
 * Perform minimization for a certain set of test values                   * 
 ***************************************************************************/
double ChiNPWrapper(glb_params base_values, double th12, double th13, double th23,
                    double delta, double dm21, double dm31, int hierarchy,
                    glb_params fit_values)
{
  double result;
  glb_params tv     = glbAllocParams();      /* Test values             */
  glb_params cv     = glbAllocParams();      /* Central values          */
  glb_params old_cv = glbAllocParams();      /* Previous central values */
  glbCopyParams(base_values, tv);
  glbGetCentralValues(cv);
  glbGetCentralValues(old_cv);

  if (!isnan(th12))
    glbSetOscParams(tv, th12, GLB_THETA_12);
  if (!isnan(th13))
    glbSetOscParams(tv, th13, GLB_THETA_13);
  if (!isnan(th23))
    glbSetOscParams(tv, th23, GLB_THETA_23);
  if (!isnan(delta))
    glbSetOscParams(tv, delta, GLB_DELTA_CP);
  if (!isnan(dm21))
    glbSetOscParams(tv, dm21, GLB_DM_21);
  if (!isnan(dm31))
    glbSetOscParams(tv, dm31, GLB_DM_31);

  if (hierarchy == HIERARCHY_NORMAL)
  {
    glbSetOscParams(tv, fabs(glbGetOscParams(tv,GLB_DM_31)), GLB_DM_31);
    glbSetOscParams(cv, fabs(glbGetOscParams(cv,GLB_DM_31)), GLB_DM_31);
    glbSetCentralValues(cv);
    result = glbChiNP(tv, fit_values, GLB_ALL);
  }
  else if (hierarchy == HIERARCHY_INVERTED)
  {
    glbSetOscParams(tv, -fabs(glbGetOscParams(tv,GLB_DM_31)), GLB_DM_31);
    glbSetOscParams(cv, -fabs(glbGetOscParams(cv,GLB_DM_31)), GLB_DM_31);
    glbSetCentralValues(cv);
    result = glbChiNP(tv, fit_values, GLB_ALL);
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
