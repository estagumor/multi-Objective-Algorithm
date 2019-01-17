/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"


/*  Test problem ZDT3
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

int nreal=30;
int nbin=0;
int ncon=0;
int nobj=2;

double
        min_realvar[30] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        max_realvar[30] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};


void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i;
    int n=30;
    double tmp=0.0;
    double g,h;
    obj[0] = xreal[0];
    for (i=1; i<n; i++)
    {
        tmp += xreal[i];
    }
    g=1+((9*tmp)/(n-1));
    h=1-sqrt(xreal[0]/g)-(xreal[0]/g)*sin(10*PI*xreal[0]);
    obj[1] = g*h;
    return;
}

