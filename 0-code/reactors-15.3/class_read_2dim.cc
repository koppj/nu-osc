#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*****************************************************
 *  class to read a 2 dim table from a data file and 
 *  perform linear interpolation on the grid
 *  the grid has to be rectangular and equally spaced
 *****************************************************/

namespace ns_reactor{

class Read_2dim{
 public:
  // at initialization provide number of x and y values and file name
  // set x_first to true if x changes and y is constant (not checked!!)
  // col gives the column of the file of the z value
  Read_2dim(const int nx, const int ny, const char *file_name, bool x_first = false);
  ~Read_2dim(void);
  double interp(const double x, const double y);

 private:
  double **data;
  double xmin, xmax, ymin, ymax, dx, dy;
};


double Read_2dim::interp(double x, double y)
{
  if(x > xmax || y > ymax || x < xmin || y < ymin){
    fprintf(stderr,"[Read_2dim] values out of range\n");
    fprintf(stderr,"x=%e, y=%e, xmax=%e, ymax=%e, xmin=%e, ymin=%e\n", 
	    x, y, xmax, ymax, xmin, ymin);
    exit(1);
  }

  if(x == xmax) x = x - 1.e-10 * fabs(x);
  if(y == ymax) y = y - 1.e-10 * fabs(y);

  int ix = int( (x - xmin) / dx );
  int iy = int( (y - ymin) / dy );

  const double fx = (x - xmin - dx * ix) / dx;
  const double fy = (y - ymin - dy * iy) / dy;

  return data[ix][iy] +
    (data[ix+1][iy] - data[ix][iy]) * fx +
    (data[ix][iy+1] - data[ix][iy]) * fy +
    (data[ix+1][iy+1] + data[ix][iy] - data[ix+1][iy] - data[ix][iy+1]) * fx * fy;
}


Read_2dim::Read_2dim(const int nx, const int ny, const char *file_name, bool x_first)
{
  /* allocate pointers to rows */  
  data = (double **) malloc( (size_t)( nx * sizeof(double*) ) );
  
  if (!data) { 
    fprintf(stderr,"[Read_2dim] allocation failure 1\n");
    exit(1);
  }

  /* allocate rows and set pointers to them */  
  data[0] = (double *) malloc((size_t)( nx * ny * sizeof(double)));
  
  if (!data[0]){
    fprintf(stderr, "[Read_2dim] allocation failure 2\n");
    exit(1);
  }

  for(int i = 1; i < nx; i++) 
    data[i] = data[i-1] + ny;

  /* read the data from file */
  FILE *fp = fopen(file_name, "r");
  if(fp == NULL){
    fprintf(stderr,"[Read_2dim] cannot open file %s\n", file_name);
    exit(1);
  }


  if(!x_first){

    for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny; j++){

        double x, y;
        if(fscanf(fp, "%lf %lf %lf\n", &x, &y, &data[i][j]) != 3){
          fprintf(stderr,"[Read_2dim] cannot read data from file %s\n", file_name);
          exit(1);
        }

        if(i == 0 && j == 0){ 
	  xmin = x; ymin = y;
	}
        if(i == nx-1 && j == ny-1){
	  xmax = x; ymax = y;
        }
      }
    }
  }else{

    for(int j = 0; j < ny; j++){
      for(int i = 0; i < nx; i++){

        double x, y;
        if(fscanf(fp, "%lf %lf %lf\n", &x, &y, &data[i][j]) != 3){
          fprintf(stderr,"[Read_2dim] cannot read data from file %s\n", file_name);
          exit(1);
        }

        if(i == 0 && j == 0){ 
	  xmin = x; ymin = y;
	}
        if(i == nx-1 && j == ny-1){
	  xmax = x; ymax = y;
        }
      }
    }
  }
  fclose(fp);

  dx = (xmax - xmin)/(nx - 1.);
  dy = (ymax - ymin)/(ny - 1.);

  if(dx <= 0. || dy <= 0.){
    fprintf(stderr,"[Read_2dim] error with boundaries (file %s)\n", file_name);
    exit(1);
  }

  return;
}

Read_2dim::~Read_2dim(void)
{
  free( (char*) (data[0]) );
  free( (char*) data );
  return;
}

}
