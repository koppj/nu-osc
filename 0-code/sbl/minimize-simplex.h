#include <cmath>
#include <cstdlib>
#include <algorithm>

#define VCT_N N_PARAM

/*************************************************************/
/*                 class vector                              */
/*************************************************************/

class Vector
{
  public:
   double x[VCT_N];
   
   // vector operators
   Vector operator + (Vector A){     
     Vector C;
     for(int i = 0; i < VCT_N; i++)
       C.x[i] = x[i] + A.x[i];
     return C;
   }      
   Vector operator - (Vector A){     
     Vector C;
     for(int i = 0; i < VCT_N; i++)
       C.x[i] = x[i] - A.x[i];
     return C;
   }   
   double operator * (Vector A){
     double p = 0.;
     for(int i = 0; i < VCT_N; i++)
       p += x[i] * A.x[i];
     return p;
   }   
   Vector operator * (double a){
     Vector C;
     for(int i = 0; i < VCT_N; i++)
       C.x[i] = a * x[i];
     return C;
   }

   // some simple vector functions   
   void set(double *a){
     for(int i = 0; i < VCT_N; i++)
       x[i] = a[i];
   } 
   void get(double *a){
     for(int i = 0; i < VCT_N; i++)
       a[i] = x[i];
   }
   void print(FILE *fp = stdout){
     for(int i = 0; i < VCT_N; i++)
       //fprintf(fp, "%8.3e\n", x[i]);
       fprintf(fp, " %e ", x[i]);
   }
   void save(char *name){
     FILE *fp = fopen(name, "w");
     if(!fp){
       fprintf(stderr, "cannot open file %s\n", name);
       exit(1);
     }
     for(int i = 0; i < VCT_N; i++)
       fprintf(fp, "%.20e\n", x[i]);
     fclose(fp);
   }
   void read(char *name){
     FILE *fp = fopen(name, "r");
     if(!fp){
       fprintf(stderr, "cannot open file %s\n", name);
       exit(1);
     }
     for(int i = 0; i < VCT_N; i++){
       if(fscanf(fp, "%lf", &x[i]) != 1){
	 fprintf(stderr, "cannot read vector from file %s\n", name);
	 exit(1);
       }
     }
     fclose(fp);
   }     
};



/*************************************************************/
/*                 class VectD                               */
/*************************************************************/

// 28/Jun/2003
// $Revision: 1.1.1.1 $
// $Date: 2005/06/01 12:03:58 $
// $Author: shindou $ 
/************************************************************************
 *  Real vector class
 *  for operation between vectors
 *  each elements of this vector is double
 *  Some functions and algorithm are from the book "Scientific C++".
 *
 *
 * ***Constructors by VectD : see compvectd.cc for more detail***
 *  Vecd V();// default
 *  Vecd V(B); // copy construction
 *  Vecd V(m);// set size of V (m) and each elements (0,0)
 *  Vecd V("filename");// input from "filename".
 *  			  // the format is as folows
 *  			  //-------------------------
 *			  // dim
 *			  // V(1)
 *			  // V(1)
 *			  // .....
 *
 * *** Member functions
 * Dim();// return dim
 * *** Operator
 * "="
 * V(i)=x;// Set V(i) to x
 * x=V(i);// return i-th element
 * V1=V2
 * "+="
 * V1+=V2: V1=V1+V2
 * "-="
 * V1-=V2: V1=V1-V2
 *  						by T. Shindou
 ************************************************************************/
/*************Caution for notations***************************************
 * vec[i] are 
 * 
 * vec[0]
 * vec[1]
 * .
 * .
 * i.e. i=0,....,dim-1
 * However vec(k)=vec[k-1]
 * So not confuse!!!
 **************************************************************************/
/**************Caution about Row- and Column- vectors**********************
 * In this class, the row- and column- vector are represented as same.
 ***************************************************************************/

class VectD
{
	private:
		int dim;// dimensions of vector
		double *vec;// elements
		
		// initialization
		inline void Init(int ndim);
		
		// deinitialisaztion
		inline void Deinit(void);

		// Preparation for copying
		inline void PrepCpy(int ndim);
	public:
		//**********constructors and destructor****************
		// default constructor
		inline VectD(void);
		// copy constructor
		inline VectD(const VectD &OtherVec);

		// generic constructors
		// for given dim (all elemens are set to 0)
		inline VectD(int ndim);
		
		// set the all elements to be val
		inline VectD(int ndim, double val);

		//destructor
		inline ~VectD(void);

		//************** access to elements***************
		// return dim
		inline int Dim(void) const
		{
			return(dim);
		}

		
		// return vec[nth-1] element
		inline double Val(int nth) const
		{
			return(vec[nth-1]);
		}

		// access vec[nth-1] element
		inline double &operator ()(int nth);
		inline double operator ()(int nth) const;

		
		// assignment operator V1=V2
		inline VectD &operator =(const VectD &rhsvec);
		// operator +=
		inline VectD &operator +=(const VectD &rhsvec);
		// operator -=
		inline VectD &operator -=(const VectD &rhsvec);
		

		//***********other functiosns*********************
		// transpose
		inline VectD T(void);
		// norm of vector
		inline double norm(void);
};// end class of VectD





/*************************************************************/
/*                 class Matrix                              */
/*************************************************************/

class Matrix{
 private:
  double **x;

 public:
  int n_row, n_col;
  Matrix(int n, int m){
    n_row = n;
    n_col = m;
    const int size = n_row * n_col;
    x = new double *[n_row];
    x[0] = new double [size];
    for(int i = 1; i < n_row; i++)
      x[i] = x[i-1] + n_col;
  }
  ~Matrix(void){
    delete x[0];
    delete x;
  }
  // sets column c to vector v
  void set_col(VectD v, int c){    
    if(c < 0 || c >= n_col)
      nrerror("Matrix::set_col: col out of range\n");    
    if(v.Dim() != n_row)
      nrerror("Matrix::set_col: cannot assign vector to matrix elements\n");
    for(int i = 0; i < n_row; i++)
      x[i][c] = v(i+1);
    return;
  }
  // returns column c as vector
  VectD get_col(int c){
    if(c < 0 || c >= n_col)
      nrerror("Matrix::get_col: col out of range\n");    
    VectD v(n_row);
    for(int i = 0; i < n_row; i++)
      v(i+1) = x[i][c];
    return v;
  }
};






#ifndef NOCHECK
# define NOCHECK
#endif

/*************************************************************/
/*                 class VectD                               */
/*************************************************************/

inline double ABS(const double x){
  return x >= 0. ? x : -x;
}

inline void VectD:: Init(int ndim)
{
	if(ndim<1)
	{
		dim=0;// initialize dim
		vec=0;// set vec to null pointer
		return;
	}
	dim=ndim;
	vec=new double [ndim];
#ifndef NOCHECK
	if(!vec)
	{
		cerr << "Error(1) in VectD::Initialization!!"<<endl;
		return;
	}
#endif // NOCHECK
}

/***************************************************************
 * Shared functions for deinitialization
 * Usage: Deinit()
 ***************************************************************/
inline void VectD:: Deinit(void)
{
	if(vec==0)
		return;
	delete vec;
	return;
}
/************** Constructors and Destructors********************/
/***************************************************************
 * Default constructor
 ***************************************************************/
inline VectD::VectD(void)
{
	Init(0);
}
/****************************************************************
 * Copy constructor
 ****************************************************************/
inline VectD::VectD(const VectD &OtherVec)
{
	Init(OtherVec.dim);
	if(dim!=0)
		memcpy(vec,OtherVec.vec,dim*sizeof(double));
}
/*****************************************************************
 * Constructor for given dim
 * Usage: VectD(dim)
 *****************************************************************/
inline VectD::VectD(int ndim)
{
	Init(ndim);
	if(dim!=0)
		memset(vec,0,dim*sizeof(double));
}
/*****************************************************************
 * Constructor for given dim and val.
 * All elements set to be val.
 * Usage: VectD(ndim,val)
 *****************************************************************/
inline VectD::VectD(int ndim, double val)
{
	Init(ndim);
	if(dim !=0)
		for(int i=0;i<dim;i++)
			vec[i]=val;
	return;
}


/*****************************************************************
 * Destructor
 *****************************************************************/
inline VectD::~VectD(void)
{
	Deinit();
}

/*****************Access to elements**********************************/
/*********************************************************************
 * return vec[nth-1] element
 *********************************************************************/
inline double &VectD::operator () (int nth)
{
#ifndef NOCHECK
	if (nth <1|| nth> dim )
	{
		cerr << "VectD: "<<dim << "<"<<nth<<endl;
		cerr << "Can't access "<< nth <<"-th element"<<endl;
	}
#endif // NOCHECK
	return(vec[nth-1]);
}
inline double VectD::operator () (int nth) const
{
#ifndef NOCHECK
	if (nth <1|| nth> dim )
	{
		cerr << "VectD: "<<dim << "<"<<nth<<endl;
		cerr << "Can't access "<< nth <<"-th element"<<endl;
	}
#endif // NOCHECK
	return(vec[nth-1]);
}

/**********************************************************************
 * Shared function for copying vector (Private)
 * Usage: PreCpy(ndim)
 **********************************************************************/
inline void VectD::PrepCpy(int ndim)
{
	if(dim != ndim)
	{
		Deinit();
		Init(ndim);
	}
	return;
}
/**********************************************************************
 * Assignment operator =
 **********************************************************************/
inline VectD &VectD::operator =(const VectD &rhsvec)
{
	PrepCpy(rhsvec.dim);
	if(dim!=0)
		memcpy(vec,rhsvec.vec,dim*sizeof(double));
	return *this;
}
/**********************************************************************
 * Assignment operator +=
 **********************************************************************/
inline VectD &VectD::operator +=(const VectD &rhsvec)
{
#ifndef NOCHECK
	if(dim!=rhsvec.Dim())
		cerr << "Operation += is unable!"<<endl;
#endif // NOCHECK
	for (int i=0;i<dim;i++)
		vec[i]+=rhsvec.Val(i+1);
	return *this;
}
/**********************************************************************
 * Assignment operator -=
 **********************************************************************/
inline VectD &VectD::operator -=(const VectD &rhsvec)
{
#ifndef NOCHECK
	if(dim!=rhsvec.Dim())
		cerr << "Operation += is unable!"<<endl;
#endif // NOCHECK
	for (int i=0;i<dim;i++)
		vec[i]-=rhsvec.Val(i+1);
	return *this;
}
/***********************************************************************
 * ******************Other member functions*****************************
 ***********************************************************************/
// transpose V^T (return V as it is)
inline VectD VectD::T(void)
{
	VectD Res(dim);
	for(int i=1;i<=dim;i++)
	{
		Res(i)=vec[i-1];
	}
	return(Res);
}


// norm of V=||V||
inline double VectD::norm(void)
{
	double Res;
	Res=0;
	for(int i=1;i<=dim;i++)
	{
		Res=Res+ABS(vec[i-1]*vec[i-1]);
	}
	return(sqrt(Res));
}

/*******************************************************************
 * ************* Operators******************************************
 * ###Numerical operators
 * V:VectD
 * d:double
 * M:CompMatD
 * -V
 * V1+V2
 * V1-V2
 * V1*V2=V1(1)V2(1)+....
 * f*V
 * V*d
 * d+V
 * V+d
 * d+V
 * d-V
 * V-d
 * V/d
 *******************************************************************/
// -V
inline VectD operator - (VectD Vec){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)=-Vec(i);
	return(Vec);
}

// Operator between Vectors
// V1+V2
inline VectD operator + (VectD Vec1, VectD Vec2){
#ifndef NOCHECK
	if(Vec1.Dim() != Vec2.Dim())
		cerr << "Operation V1+V2 is unable! "<<endl;
#endif // NOCHECK
	for(int i=1;i<=Vec1.Dim();i++)
		Vec1(i)+=Vec2(i);
	return(Vec1);
}

// V1-V2
inline VectD operator - (VectD Vec1, VectD Vec2){
#ifndef NOCHECK
	if(Vec1.Dim() != Vec2.Dim())
		cerr << "Operation V1-V2 is unable! "<<endl;
#endif // NOCHECK
	for(int i=1;i<=Vec1.Dim();i++)
		Vec1(i)-=Vec2(i);
	return(Vec1);
}

// V1*V2
// Caution: V1*V2=V1^T.V2=V1(1)*V2(1)+V1(2)*V2(2)+....( NEQ V^dag.V2)
inline double operator * (VectD Vec1, VectD Vec2){
#ifndef NOCHECK
	if(Vec1.Dim() != Vec2.Dim())
	    cerr << "Operation V1*V2 is unable! "<<endl;
#endif // NOCHECK
	double Ans;
	Ans=0;
	for(int i=1;i<=Vec1.Dim();i++)
	{
		Ans=Ans+Vec1(i)*Vec2(i);
	}
	return(Ans);
}

// Scalar * Vector
// d*V
inline VectD operator * (double numD, VectD Vec){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)*=numD;
	return(Vec);
}
// V*d
inline VectD operator * (VectD Vec,double numD){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)*=numD;
	return(Vec);
}

// Scalar +- Vector & Vector +- Scalar
// d+V
inline VectD operator + (double numD,VectD Vec){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)+=numD;
	return(Vec);
}
// V+d
inline VectD operator + (VectD Vec,double numD){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)+=numD;
	return(Vec);
}

// d-V
inline VectD operator - (double numD,VectD Vec){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)-=numD;
	return(-Vec);
}
// V-d
inline VectD operator - (VectD Vec,double numD){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)-=numD;
	return(Vec);

}

// Vector/Scalar
// V/d
inline VectD operator / (VectD Vec, double numD){
	for(int i=1;i<=Vec.Dim();i++)
		Vec(i)/=numD;
	return(Vec);
}

/****************************************************************
 *************** Another Functions********************************
 *
 * transpose(V);// transpose (return V as it is)
 * norm(V); // return norm of V
 * normalize(V);// return normalized vector 
 ****************************************************************/
inline VectD transpose(VectD Vec){
	return(Vec);
}
inline double norm(VectD Vec){
	return(sqrt(ABS(Vec*Vec)));
}
inline VectD normalize(VectD Vec){
	return(Vec/norm(Vec));
}
/****************I/O stream***************************************
 * The data format is 
 *
 * ndim
 * v(1) v(2) v(3) .....
 *****************************************************************
inline ostream &operator <<(ostream &ost,const VectD &outvec)
{
	ost << outvec.Dim()<<endl;
	for(int i=1;i<=outvec.Dim();i++)
		ost << outvec.Val(i)<< ' ';
	ost << endl;
	return(ost);
}
inline istream &operator >>(istream &ist, VectD &invec)
{
	int ndim;
	ist >> ndim;
	if(ist.bad()) return(ist);// Check the iostream
	VectD tmp(ndim);// Temporary
	for(int i=1;i<=tmp.Dim();i++){
		ist >> tmp(i);
		if(ist.bad()) return(ist);// Check the iostream
	}
#ifndef NOCHECK
	if(invec.Dim()!=tmp.Dim()){
		cerr << "Istream: Dimesnion of vector is incorrect!!"<<endl;
		exit(1);
	}
#endif // NOCHECK
	invec=tmp;
	return(ist);
}
*/

