/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"


#define MYSIGN(x) ((x)>0?1.0:-1.0)

/*  Test problem CEC09 CF6
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

int nreal=4;
int nbin=0;
int ncon=2;
int nobj=2;

double
        min_realvar[4] = {0.0,-2.0,-2.0,-2.0},
        max_realvar[4] = {1.0,2.0,2.0,2.0};


void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
                unsigned int j;
                double sum1, sum2, yj;

                sum1   = sum2   = 0.0;
                for(j = 2; j <= nreal; j++)
                {
                        if (j % 2 == 1)
                        {
                                yj     = xreal[j-1] - 0.8*xreal[0]*cos(6.0*PI*xreal[0] + j*PI/nreal);
                                sum1  += yj*yj;
                        }
                        else
                        {
                                yj     = xreal[j-1] - 0.8*xreal[0]*sin(6.0*PI*xreal[0] + j*PI/nreal);
                                sum2  += yj*yj;
                        }
                }
                obj[0] = xreal[0] + sum1;
                obj[1] = (1.0 - xreal[0])*(1.0 - xreal[0]) + sum2;
                constr[0] = xreal[1]-0.8*xreal[0]*sin(6.0*xreal[0]*PI+2.0*PI/nreal) - MYSIGN((xreal[0]-0.5)*(1.0-xreal[0]))*sqrt(fabs((xreal[0]-0.5)*(1.0-xreal[0])));
                constr[1] = xreal[3]-0.8*xreal[0]*sin(6.0*xreal[0]*PI+4.0*PI/nreal) - MYSIGN(0.25*sqrt(1-xreal[0])-0.5*(1.0-xreal[0]))*sqrt(fabs(0.25*sqrt(1-xreal[0])-0.5*(1.0-xreal[0])));

    return;
}

