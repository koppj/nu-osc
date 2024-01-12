#include "class_minimize.h"

double Minimize::find_min(double *x, bool *fixed)
{
   double w = minimize(x, fixed);
   double v = 2.*w;

   // restart until min converges
   while(fabs(v - w) > fabs(w) * 0.01 + 1.e-7){
      v = w;
      w = minimize(x, fixed);
      //if(w > v) 
      //fprintf(stderr, "WARNING [minimize]: new minimum (%e) larger than previous one (%e)\n", w, v);
   }
   if(fabs(w) < 1.e-5)
     fprintf(stderr, "WARNING [minimize]: minimum smaller than 1.e-5\n");
   return w;
}

double Minimize::find_min(double *x)
{
   bool fixed[N_PAR_MIN];
   for(int i = 0; i < N_PAR_MIN; i++)
     fixed[i] = false;
   return find_min(x, fixed);
}



/************************************************* 
 * 1-dim line minimization along the direction n
 *************************************************/

double Minimize::line_min(double *P, double *n, int *num)
{   
  double P0[N_PAR_MIN];
  for(int i = 0; i < N_PAR_MIN; i++) 
    P0[i] = P[i];

  double w1 = min_func(P);
  double lambda = 1.;

  // try a first step
  for(int i = 0; i < N_PAR_MIN; i++) 
    P[i] = P0[i] + n[i] * lambda;
  double w2 = min_func(P);

  int m = 2;

  // check if step size is much too big
  while(w2 - w1 > 5.*w1){
    lambda *= 0.5;	
    for(int i = 0; i < N_PAR_MIN; i++) 
      P[i] = P0[i] + n[i] * lambda;
    w2 = min_func(P);
    m++;
  }
        
  if(w2 > w1){      
    lambda *= -1;   // change direction
    for(int i = 0; i < N_PAR_MIN; i++) 
      P[i] = P0[i];
    w2 = w1;
  }      
  // w2 is the smallest point found so far

  // search for the region of the minimum
  do{
    w1 = w2;
    
    for(int i = 0; i < N_PAR_MIN; i++){ 
      P[i] += n[i] * lambda;    
      if(fabs(P[i]) > 1.e30)
	fprintf(stderr, "MINIMIZE: WARNING: parameter gets very big %e %d!\n", P[i], i);
    }
    w2 = min_func(P);

    if(++m == MAX_N_LINE_MIN)
      error("[line_min]: no minimum found");
    
    if(w1 - w2 > w2)
      lambda *= 0.5;
    else{
      if(w1 - w2 < 0.1 * w2 && w1 -w2 > 0.)
        lambda *= 1.5;
    }

  }while(w1 - w2 > 0.);

  // w1 is the smallest point found so far at
  for(int i = 0; i < N_PAR_MIN; i++) 
    P0[i] = P[i] - n[i] * lambda;

  // starting the loop around the minimum
  lambda *= .5;
  double delta;
  do{
    for(int i = 0; i < N_PAR_MIN; i++) 
      P[i] = P0[i] + n[i] * lambda;
    w2 = min_func(P);
    delta = w2 - w1;

    if(++m == MAX_N_LINE_MIN)
      fprintf(stderr, "too many steps in line_min");

    if(w2 > w1) 
      lambda *= -2./3.;
    else{
      for(int i = 0; i < N_PAR_MIN; i++) 
        P0[i] = P[i];
      w1 = w2;
    }

  }while(fabs(delta) > ACC * (fmin(w2,w1) + .1));

  //fprintf(stderr, "%d calls to chisq\n", m);
  *num += m;

  return fmin(w2, w1);
}

/*******************************************************
 * minimization by Powells algorithm
 * modified according to Numerical Recipes chapt. 10.5
 *
 * *P contains the starting point and 
 * after returning *P contains the minimum point  
 *******************************************************/

double Minimize::minimize(double *P, bool *fixed)
{
  // initializing the directions
  double n[N_PAR_MIN][N_PAR_MIN];
   
  for(int i = 0; i < N_PAR_MIN; i++){
    for(int j = 0; j < N_PAR_MIN; j++)
      n[i][j] = 0.;

    if(P[i] == 0.) 
      n[i][i] = 0.1;
    else
      n[i][i] = P[i] * .2;
  }
   
  int first, last;
  for(first = 0; fixed[first]; first++);
  for(last = N_PAR_MIN-1; fixed[last]; last--);
  if(first > last)
    error("[minimize] first > last\n");
  
  int chisq_calls = 0;
  if(first == last){
    // only one-dimensional minimization
    return line_min(P, n[first], &chisq_calls);
  }

  double Q[N_PAR_MIN][N_PAR_MIN];
  
  for(int i = 0; i < N_PAR_MIN; i++)
    Q[last][i] = P[i];
    
  double f2 = min_func(P);
  chisq_calls = 1;

  int mm = 0;
  for(int m = 0; m < MAX_N_MIN; m++, mm++){

    // the initial point
    double P0[N_PAR_MIN];
    
    for(int i = 0; i < N_PAR_MIN; i++)
      P0[i] = Q[last][i];
    
    double f1, f0 = f2; 

    // the first line minimization
    f2 = line_min(P0, n[first], &chisq_calls);
    
    for(int i = 0; i < N_PAR_MIN; i++)
      Q[first][i] = P0[i];
    
    double delta_max = f0 - f2; 
    int i_max = 0;

    // all the other line minimizations
    for(int i = first + 1; i <= last; i++){
      if(!fixed[i]){

	int previous;
        for(previous = i-1; fixed[previous]; previous--);

        f1 = f2;
        f2 = line_min(Q[previous], n[i], &chisq_calls);
	
	for(int j = 0; j < N_PAR_MIN; j++)
          Q[i][j] = Q[previous][j];

        if(f1 - f2 > delta_max){
	  delta_max = f1 - f2;
	  i_max = i;
	}
      }
    }

    //    fprintf(stderr, "%e\n", f2);

    // minimum found?
    if(fabs(f2 - f0) < ACC * (f2 + 0.1) ){
      for(int i = 0; i < N_PAR_MIN; i++)
        P[i] = Q[last][i];
      /*
      fprintf(stderr, "%d loops in minimize -> %d calls to line_min, ", 
              m, m * (N_PAR_MIN+1) - 1);
      fprintf(stderr, "%d calls to chisq\n", chisq_calls);
       */         
      return f2;
    }
  }
  fprintf(stderr, "too many loops in minimize");  
  fprintf(stderr, "%f\n", f2);
  return f2;
}


   

