#include <cmath>

#define EIGEN_EPSILON_VAL    1.0e-15
#define EIGEN_EPSILON_NORM   0.
#define EIGEN_EPSILON_VALQ  (EIGEN_EPSILON_VAL * EIGEN_EPSILON_VAL)
#define EIGEN_EPSILON_NORMQ (EIGEN_EPSILON_NORM * EIGEN_EPSILON_NORM)

#define PRECISION 1e-11

inline double norm(double x)
{
   return x*x;
}

/************************************************/
/* Real Symmetric matrix: (Nt)(B)(N) -> diag(L) */
/************************************************/

bool eigenvectors(const int size, void *B_in, void *N_out, void *L_out)
{
   double A[size][size];
   double (&B)[size][size] = *(double (*)[size][size]) B_in;
   double (&N)[size][size] = *(double (*)[size][size]) N_out;
   double (&L)[size] = *(double (*)[size]) L_out;

   double nrmx = 0., notd = 0., subtract;
   register double tmp, tmq;

   /* Evaluate the norm of the matrix */

   for(int p = 0; p < size; p++) {

      N[p][p] = 1.;
      nrmx += norm(A[p][p] = B[p][p]);

      for(int q = p+1; q < size; q++) {

	 N[p][q] = N[q][p] = 0.;
	 notd += 2. * norm(A[p][q] = B[p][q]);
      }
   }

   nrmx += notd;

   /* Diagonalize matrix and find eigenvectors */

   do {
      subtract = 0.;

      for(int p = 0; p < size-1; p++)
	for(int q = p+1; q < size; q++) {

	   double scale_n = 2. * norm(A[p][q]);
	   double scale_t = EIGEN_EPSILON_VALQ * (scale_n + norm(A[p][p]) + norm(A[q][q]));

	   if(scale_n <= scale_t)
	     continue;

	   subtract += scale_n;

	   tmp = atan2(2.* A[p][q], A[q][q] - A[p][p]) / 2.;
	   double cos_n = cos(tmp);
	   double sin_n = sin(tmp);

	   tmp = cos_n*cos_n*A[p][p] + sin_n*sin_n*A[q][q] - 2.*cos_n*sin_n*A[p][q];
	   tmq = sin_n*sin_n*A[p][p] + cos_n*cos_n*A[q][q] + 2.*cos_n*sin_n*A[p][q];
	   A[p][q] = A[q][p] = 0.;
	   A[p][p] = tmp;
	   A[q][q] = tmq;

	   for(int r = 0; r < size; r++) {

	      tmp = N[r][p];
	      tmq = N[r][q];
	      N[r][p] = cos_n * tmp - sin_n * tmq;
	      N[r][q] = sin_n * tmp + cos_n * tmq;

	      if(r < p) {

		 tmp = A[r][p];
		 tmq = A[r][q];
		 A[r][p] = cos_n * tmp - sin_n * tmq;
		 A[r][q] = sin_n * tmp + cos_n * tmq;
	      } 
	      else if(p < r && r < q) {

		 tmp = A[p][r];
		 tmq = A[r][q];
		 A[p][r] = cos_n * tmp - sin_n * tmq;
		 A[r][q] = sin_n * tmp + cos_n * tmq;
	      }
	      else if(q < r) {

		 tmp = A[p][r];
		 tmq = A[q][r];
		 A[p][r] = cos_n * tmp - sin_n * tmq;
		 A[q][r] = sin_n * tmp + cos_n * tmq;
	      }
	   }
	}

      notd -= subtract;

   } while(subtract && fabs(notd/nrmx) >= EIGEN_EPSILON_NORMQ);

   for(int p = 0; p < size; p++)
     L[p] = A[p][p];

#ifdef DEBUG

   /* Check whether V is really an orthogonal matrix */

   double err1 = 0.;

   for(int p = 0; p < size; p++)
     for(int q = 0; q < size; q++) {

	tmp = 0.;
	for(int r = 0; r < size; r++)
	  tmp += N[r][p] * N[r][q];

	err1 += norm(tmp - double(p==q));
     }

   /* Check whether (N+)(B)(N) is really diag(L) */

   double err2 = 0.;

   for(int p = 0; p < size; p++)
     for(int q = 0; q < size; q++) {

	tmp = 0.;
	for(int r = 0; r < size; r++)
	  tmp += N[p][r] * L[r] * N[q][r];

	err2 += norm(B[p][q] - tmp);
     }

   if(sqrt(err1/size + err2/nrmx) >= PRECISION)
     return true;

#endif

   return false;
}

/********************************/
/* Invert a matrix: Inv[B] -> C */
/********************************/

bool invert(const int size, void *B_in, void *C_out)
{
   double N[size][size], L[size];
   double (&B)[size][size] = *(double (*)[size][size]) B_in;
   double (&C)[size][size] = *(double (*)[size][size]) C_out;

   if(eigenvectors(size, B, N, L))
     return true;

   for(int p = 0; p < size; p++) {
      if(L[p] == 0.) return true;
      L[p] = 1. / L[p];
   }

   for(int p = 0; p < size; p++)
     for(int q = 0; q < size; q++) {
	C[p][q] = 0.;
	for(int r = 0; r < size; r++)
	  C[p][q] += N[p][r] * N[q][r] * L[r];
     }

#ifdef DEBUG

   /* Check whether C is really the inverse of B */

   double err = 0.;

   for(int p = 0; p < size; p++)
     for(int q = 0; q < size; q++) {

	double tmp = 0.;
	for(int r = 0; r < size; r++)
	  tmp += B[p][r] * C[r][q];

	err += norm(tmp - double(p==q));
     }

   if(sqrt(err/size) >= PRECISION)
     return true;

#endif

   return false;
}

/********************************************/
/* Solve any solvable linear system: B -> L */
/********************************************/

bool singsolve(const int size, void *B_in, void *L_out)
{
   double A[size][size+1];
   const int size1 = size + 1;
   double (&B)[size][size1] = *(double (*)[size][size1]) B_in;
   double (&L)[size] = *(double (*)[size]) L_out;
   int sort[size];

   double nrmx = 0.;
   register double tmp;

   /* Evaluate the norm of the matrix */

   for(int r = 0; r < size; r++) {
      sort[r] = r;

      for(int c = 0; c <= size; c++)
	nrmx += norm(A[r][c] = B[r][c]);
   }

   /* Step 1: Row-reduce the matrix */

   int pv;
   for(pv = 0; pv < size; pv++) {

      /* Find the largest element in the matrix */

      tmp = 0.;
      int nr = pv;
      int nc = pv;

      for(int r = pv; r < size; r++)
	for(int c = pv; c < size; c++)
	  if(fabs(A[r][c]) > tmp) {
	     tmp = fabs(A[r][c]);
	     nr = r;
	     nc = c;
	  }

      /* Check whether the matrix is singular */
      
      if(tmp < PRECISION)
	break;

      /* Exchange the current row with the row containing the largest element */

      if(pv < nr)
	for(int c = pv; c <= size; c++) {
	   tmp = A[pv][c];
	   A[pv][c] = A[nr][c];
	   A[nr][c] = tmp;
	}

      /* Exchange the current column with the column containing the largest element */

      if(pv < nc) {
	 for(int r = 0; r < size; r++) {
	    tmp = A[r][pv];
	    A[r][pv] = A[r][nc];
	    A[r][nc] = tmp;
	 }
	 int p = sort[pv];
	 sort[pv] = sort[nc];
	 sort[nc] = p;
      }

      /* Execute Gauss reduction */

      for(int r = pv+1; r < size; r++) {
	 tmp = A[r][pv] / A[pv][pv];
	 for(int c = pv; c <= size; c++)
	   A[r][c] -= A[pv][c] * tmp;
      }
   }

   /* Step 2: Solve the linear system */

   for(int p = pv; p < size; p++)
     L[p] = 0.;

   for(pv--; pv >= 0; pv--) {

      /* Evaluate the p-th component of the solution */

      tmp = L[pv] = -A[pv][size] / A[pv][pv];
 
      /* Subtract the p-th component of the solution from the second term */

      for(int r = 0; r < pv; r++)
	A[r][size] += A[r][pv] * tmp;
   }

   /* Step 3: Sort the components of the solution */

   for(int p = 0; p < size; p++) {

      int q;
      while(q = sort[p], q != p) {

	 sort[p] = sort[q];
	 sort[q] = q;

	 tmp = L[p];
	 L[p] = L[q];
	 L[q] = tmp;
      }
   }

   /* Step 4: Check the result and return the error */

   double err = 0.;
   for(int r = 0; r < size; r++) {

      tmp = B[r][size];
      for(int c = 0; c < size; c++)
	tmp += B[r][c] * L[c];

      err += norm(tmp);
   }

   return bool(sqrt(err/nrmx) >= PRECISION);
}
