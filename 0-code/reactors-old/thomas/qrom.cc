/**********************************************************/
/*            integration routines                        */
/**********************************************************/

#include "qrom.h"

void nrerror(const char error_text[])
{
   fprintf(stderr, "%s\n", error_text);
   exit(1);
}


#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n, int k)
{
	double x,tnm,sum,del;
	static double s[3];
	int it,j;

	if (n == 1) {
		return (s[k]=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s[k]=0.5*(s[k]+(b-a)*sum/tnm);
        
		return s[k];
	}
}

#define NR_END 1
#define FREE_ARG char*

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{  
   void nrerror(const char error_text[]);
   double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
        void nrerror(const char error_text[]);
   
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
          for (i=1;i<=n-m;i++) {
	     ho=xa[i]-x;
	     hp=xa[i+m]-x;
	     w=c[i+1]-d[i];
	     if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
	     den=w/den;
	     d[i]=hp*den;
	     c[i]=ho*den;
	  }
	  *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}




#define JMAX 20
#define JMAXP (JMAX+1)
#define K 4


double qromb1(double (*func)(double), double a, double b, double eps, int k)
{
    void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
    double trapzd(double (*func)(double), double a, double b, int n, int k);
    void nrerror(const char error_text[]);
    double ss,dss;
    double s[JMAXP+1],h[JMAXP+1];
    int j;

  
    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
      s[j]= trapzd(func,a,b,j,k);
      if (j >= K) {
	polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
	
	if(fabs(dss) <= eps*fabs(ss)) 
	  return ss;
	if(fabs(dss)<eps && fabs(ss)<10.0*eps){
	  //fprintf(stderr,"qromb1 WARNING: absolute accuracy: %e!\n", ss);
	  return ss;
	}
      }
      s[j+1]=s[j];
      h[j+1]=0.25*h[j];
    }
    nrerror("Too many steps in routine qromb1");
    return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K

/************************************************************************
 * routines for integrating without evaluating at boundaries
 ************************************************************************/

double midpnt(double (*func)(double), double a, double b, int n, int k)
{
  static double s[3];
  int it, j;

  if(k >= 3) nrerror("midpnt: k out of range"); 

  if(n == 1)
    return (s[k] = (b - a) * FUNC(0.5 * (a + b)));

  else{
    for(it = 1, j = 1; j <= n-1; j++) it *= 3;

    double tnm = it;
    double del = (b - a) / (3. * tnm);
    double ddel = del + del;
    double x = a + 0.5 * del;
    double sum = 0.;
    for(j = 1; j <= it; j++){

      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    s[k] = (s[k] + (b - a) * sum / tnm)/3.;
    return s[k];
  }
}



#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

double qromo(double (*func)(double), double a, double b, double eps, int k)
{
    void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
    double midpnt(double (*func)(double), double a, double b, int n, int k);
    void nrerror(const char error_text[]);
    double ss,dss;
    double s[JMAXP+1],h[JMAXP+1];
    int j;

  
    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
      s[j]= midpnt(func,a,b,j,k);
      if (j >= K) {
	polint(&h[j-K],&s[j-K],K,0.,&ss,&dss);
	
	if(fabs(dss) <= eps*fabs(ss)) 
	  return ss;
	if(fabs(dss)<eps && fabs(ss)<10.0*eps){
	  //fprintf(stderr,"qromb1 WARNING: absolute accuracy: %e!\n", ss);
	  return ss;
	}
      }
      s[j+1]=s[j];
      h[j+1]=h[j]/9.;
    }
    nrerror("Too many steps in routine qromo");
    return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K

/************************************************************************
 * routines for integrating from a > 0 to infinity
 ************************************************************************/

// using: \int_a^\infty dx f(x) = \int_0^{1/a} dt 1./t^2 f(1/t) 

double (*func_glob0)(double);
double func_transf0(double t){
  return (*func_glob0)(1./t) / (t*t);
}
double (*func_glob1)(double);
double func_transf1(double t){
  return (*func_glob1)(1./t) / (t*t);
}
double (*func_glob2)(double);
double func_transf2(double t){
  return (*func_glob2)(1./t) / (t*t);
}

double qromo_inf(double (*func)(double), double a, double eps, int k)
{
  if(k == 0){
    func_glob0 = func;
    return qromo(func_transf0, 0., 1./a, eps, k);
  }else if(k == 1){
    func_glob1 = func;
    return qromo(func_transf1, 0., 1./a, eps, k);
  }else{
    func_glob2 = func;
    return qromo(func_transf2, 0., 1./a, eps, k);
  }
}
