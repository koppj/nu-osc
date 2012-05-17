#include "definitions.h"

#define N_PARAM 7

#include "minimize-simplex.h"

#define ACC 1.e-5
#define MAX_N_MIN 200000


/*******************************************************
 * minimization by the downhill-simplex algorith
 * see Numerical Recipes chapt. 10.4
 *
 * *P contains the starting point and 
 * after returning *P contains the minimum point  
 *******************************************************/

// global variables
bool fixed_global[N_PARAM], incl_gl[NUM_EXP];
int n_min;

// indices for the parameters
#define Ue4 0
#define Ue5 1
#define Um4 2
#define Um5 3
#define DEL 4
#define D41 5
#define D51 6

// converting from params struct to Vector
void param2vec(params *p, Vector *P){
  P->x[Ue4] = log(p->Ue[I4]);
  P->x[Ue5] = log(p->Ue[I5]);
  P->x[Um4] = log(p->Um[I4]);
  P->x[Um5] = log(p->Um[I5]);
  P->x[DEL] = p->delta;
  P->x[D41] = log(p->dmq[I4]);
  P->x[D51] = log(p->dmq[I5]);
  return;
}
void vec2param(Vector *P, params *p){
  p->Ue[I4]  =  exp(P->x[Ue4]);
  p->Ue[I5]  =  exp(P->x[Ue5]);
  p->Um[I4]  =  exp(P->x[Um4]);
  p->Um[I5]  =  exp(P->x[Um5]);
  p->delta   =  P->x[DEL];
  p->dmq[I4] =  exp(P->x[D41]);
  p->dmq[I5] =  exp(P->x[D51]);
  return;
}

double minimize_simplex(Vector *P, int *n_evaluation);
double minimize_simplex(Vector *P, const bool fixed[N_PARAM]);

// the minimization routine visible to main
double min_chisq(params *p, const fix_params f, const bool incl[NUM_EXP])
{
  Vector P;
  param2vec(p, &P);

  bool fixed[N_PARAM];
  fixed[Ue4] = f.Ue[I4];
  fixed[Ue5] = f.Ue[I5];
  fixed[Um4] = f.Um[I4];
  fixed[Um5] = f.Um[I5];
  fixed[DEL] = f.delta;
  fixed[D41] = f.dmq[I4];
  fixed[D51] = f.dmq[I5];

  for(int i = 0; i < NUM_EXP; i++)
    incl_gl[i] = incl[i];

  const double w = minimize_simplex(&P, fixed);
  vec2param(&P, p);

  while(p->delta > 2*M_PI)
    p->delta -= 2*M_PI;
  while(p->delta < 0.)
    p->delta += 2*M_PI;

  return w;
}

// the chisq called by the minimisor
#define BIG 1.e6

double chisq(Vector P){
  params p;
  vec2param(&P, &p);

  // check boundaries
  const double de4 = sqr(p.Ue[I4]);
  const double de5 = sqr(p.Ue[I5]);
  const double dm4 = sqr(p.Um[I4]);
  const double dm5 = sqr(p.Um[I5]);

  if(de4 + de5 > 1.)
    return BIG * (1. + sqr(de4 + de5 - 1.));
  if(dm4 + dm5 > 1.)
    return BIG * (1. + sqr(dm4 + dm5 - 1.));
  if(de4 + dm4 > 1.)
    return BIG * (1. + sqr(de4 + dm4 - 1.));
  if(de5 + dm5 > 1.)
    return BIG * (1. + sqr(de5 + dm5 - 1.));

  const double dmqmax = pow(10., LDM2MAX);

  if(p.dmq[I5] > dmqmax)
    return BIG * (1. + sqr(p.dmq[I5] - dmqmax));

#if defined(CHISQ) && !defined(LOOP_DMQ)  
  if(p.dmq[I4] > p.dmq[I5])
    return BIG * (1. + sqr(p.dmq[I4] - p.dmq[I5]));  
#else
  if(p.dmq[I4] > dmqmax)
    return BIG * (1. + sqr(p.dmq[I4] - dmqmax));  
#endif
    
  return chisq_main(p, incl_gl);
}


// translating from Vector(N_PARAM) to VectD(n_min)
VectD copy_vect(const Vector P)
{
  VectD v(n_min);
  int j = 0;  
  for(int i = 0; i < N_PARAM; i++){
    if(!fixed_global[i]){
      v(j+1) = P.x[i];
      j++;
    }
  }
  return v;
} 
Vector copy_vect(const VectD v, Vector *P)
{
  int j = 0;  
  for(int i = 0; i < N_PARAM; i++){
    if(!fixed_global[i]){
      P->x[i] = v(j+1);
      j++;
    }
  }
  return *P;
} 


/* extrapolates by a factor fact through the face of the simplex
 * across from the highest point, tries it, and replaces the high
 * point if the new point is better
 */
double amoeba(Matrix *p, VectD *y, VectD *psum, Vector *P, 
              const int i_highest, const double fact)
{
  const double fact1 = (1. - fact)/n_min;
  const double fact2 = fact1 - fact;

  VectD ptry(*psum * fact1 - p->get_col(i_highest-1) * fact2);
  const double ytry = chisq(copy_vect(ptry, P));

  if(ytry < (*y)(i_highest)){
    (*y)(i_highest) = ytry;
    *psum += ptry - p->get_col(i_highest-1);
    p->set_col(ptry, i_highest-1);
  }
  return ytry;
}
  
double minimize_simplex(Vector *P, int *n_evaluation)
{
  // initializing the directions
  Vector n[N_PARAM];
  for(int i = 0; i < N_PARAM; i++)
    for(int j = 0; j < N_PARAM; j++)
      n[i].x[j] = 0.;

  // Ue4, Ue5, Um4, Um5
  for(int i = 0; i < 4; i++)
    n[i].x[i] = 1.;

  n[DEL].x[DEL] = 0.1 * M_PI;

  n[D41].x[D41] = .3;
  n[D51].x[D51] = .3;

  // the n_min+1 points of the simplex
  Matrix p(n_min, n_min+1);
  p.set_col(copy_vect(*P), 0);
  int j = 1;  
  for(int i = 0; i < N_PARAM; i++){
    if(!fixed_global[i]){
      p.set_col(p.get_col(0) + copy_vect(n[i]), j);
      j++;
    }
  }

  // the chisq values at the n_min+1 points
  VectD y(n_min+1);
  for(int i = 0; i < n_min+1; i++)
    y(i+1) = chisq(copy_vect(p.get_col(i), P));

  *n_evaluation += n_min+1;

  // the sum of all n_min+1 vectors
  VectD psum(p.get_col(0));
  for(int i = 1; i < n_min+1; i++)
    psum += p.get_col(i);
  
  // the loop
  for(;;){

    // find the highest, next-to-highest, and lowest points
    int i_lowest = 1, i_nthighest;
    int i_highest = y(1) > y(2) ? (i_nthighest = 2, 1) : (i_nthighest = 1, 2);
    for(int i = 1; i <= n_min+1; i++){
      if(y(i) <= y(i_lowest)) 
	i_lowest = i;
      if(y(i) > y(i_highest)){
	i_nthighest = i_highest; 
	i_highest = i;
      }else if(y(i) > y(i_nthighest) && i != i_highest)
	i_nthighest = i;
    }
    
    //fprintf(stderr, "%g\n", y(i_lowest));
    
    // check if done
    if( fabs(y(i_highest) - y(i_lowest))/(fabs(y(i_lowest)) + 10. * ACC) < ACC){
      copy_vect(p.get_col(i_lowest-1), P);
      return y(i_lowest);
    }
    if(*n_evaluation >= MAX_N_MIN){
      fprintf(stderr, "WARNING: too many steps in minimize-simplex\n");
      fprintf(stderr, "  stopping at epsilon = %g\n", 
              fabs(y(i_highest) - y(i_lowest))/(fabs(y(i_lowest)) + 10. * ACC));
      copy_vect(p.get_col(i_lowest-1), P);
      return y(i_lowest);
    }

    /* begin new iteration */
   
    // first reflect the simplex from the high point
    double ytry = amoeba(&p, &y, &psum, P, i_highest, -1.);
    (*n_evaluation)++;

    if(ytry <= y(i_lowest)){
      // additional extrapolation by factor 2
      ytry = amoeba(&p, &y, &psum, P, i_highest, 2.);
      (*n_evaluation)++;
  
    }else if(ytry >= y(i_nthighest)){
      const double ysave = y(i_highest);
      // try intermediate point
      ytry = amoeba(&p, &y, &psum, P, i_highest, .5);
      (*n_evaluation)++;

      if(ytry >= ysave){
	//contract around the lowest point
	for(int i = 0; i < n_min+1; i++){

	  if(i+1 != i_lowest){
	    p.set_col(0.5 * (p.get_col(i_lowest-1) + p.get_col(i)), i);
            y(i+1) = chisq(copy_vect(p.get_col(i), P));
	  }
	  if(i == 0)
	    psum = p.get_col(i);
          else
	    psum += p.get_col(i);
	}
        *n_evaluation += n_min;
      }
    }
  }
  return 0.;
}


/****************************************************** 
 * repeating the simplex algorithm until it stabilizes 
 ******************************************************/

double minimize_simplex(Vector *P, const bool fixed[N_PARAM])
{
  // find the dimension of the minimization 
  n_min = 0;
  for(int i = 0; i < N_PARAM; i++){
    fixed_global[i] = fixed[i];
    if(!fixed[i]) n_min++;
  }

  int n_evaluate = 0;
  
  double q_old, q_new = minimize_simplex(P, &n_evaluate);
  int n_simplex = 1;
  
  do{
    q_old = q_new;
    q_new = minimize_simplex(P, &n_evaluate);
    n_simplex++;
    //fprintf(stderr, "%g  %d\n", q_old - q_new, n_evaluate);
  }while(q_old - q_new > 10. * ACC && n_evaluate < MAX_N_MIN);
  /*    
  fprintf(stderr, "%d chisq evaluations, %d calls to simplex\n", 
	  n_evaluate, n_simplex);
  */
  return q_new;
}


