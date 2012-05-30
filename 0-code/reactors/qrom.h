/**********************************************************/
/*            integration routines                        */
/**********************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace ns_reactor
{

double qromb1(double (*func)(double), double a, double b, double eps, int k=0);

/************************************************************************
 * integrating without evaluating at boundaries
 ************************************************************************/

double qromo(double (*func)(double), double a, double b, double eps, int k=0);

/************************************************************************
 * integrating from a > 0 to infinity
 ************************************************************************/

double qromo_inf(double (*func)(double), double a, double eps, int k=0);

} // end namespace

