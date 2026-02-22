/**
Moving least squares (MLS) implementation in Basilisk (2D only). */

#include "utils.h"

/**
Cubic spline weighting function. r = |x - xi| and h is the support radius.*/
double cubic_spline (double r, double h)
{
    double s = r/h;
    if (s >= 0 && s <= 0.5)
        return 2./3. - 4*sq(s) + 4*cube(s);
    else if (s > 0.5 && s <= 1)
        return 4./3. - 4*s + 4*sq(s) - 4*cube(s)/3.;
    else
        return 0;
}

/**
Gaussian weighting function*/
double gaussian_function (double r, double h)
{
    return exp(-sq(r/h));
}

void fill_basis (int m, double v[m], coord p)
{
    v[0] = 1; v[1] = p.x; v[2] = p.y;

    if (m == 6) {
        v[3] = p.x*p.y; v[4] = sq(p.x); v[5] = sq(p.y);
    }
}

int gauss_elim(int m, int n, double matrix[m][n], double sol[m]);

/**
This function finds a value at a given point pint by using MLS interpolation given np
number of points, where p and u stores the coordinates and value for each point, resp. 

If linear = true, then 1st order interpolation is used, which requires atleast 3 points.
If linear = fals, then we use quadratic interpolation, which is second order, but requires
at least 6 points.

Range is used as input into the weighting function and acts as the range of influence points
have on the point we're interpolating to. It should be based on Delta (the cell size).

One can also choose to use a cubic spline or gaussian function for weighting.*/

double interpolate_mls (int np, coord p[np], double u[np], coord pint, double range, double h,
                        bool linear = true, double (*weighting_func) (double,double) = cubic_spline)
{
    if (linear && np < 3) {
        fprintf(stderr, "WARNING: linear moving least squares interpolation requires at least 3 points"
                         ", but only %d were given\n", np);
        return nodata;
    }
    else if (!linear && np < 6) {
         fprintf(stderr, "WARNING: quadratic moving least squares interpolation requires at least 6 points"
                         ", but only %d were given\n", np);
        return nodata; 
    }

    const int m = linear? 3: 6;
    
    /** 
    Allocate memory for vectors and matrices.*/
    double (*A)[m+1] = malloc(m * sizeof(*A));
    double * w = malloc(np * sizeof(double)); // weights
    double * c = malloc(m * sizeof(double)); // coefficients

    for (int i = 0; i < m; i++) {
        //A[i] = malloc(m + 1 * sizeof(double)); // + 1 for augmented matrix
        for (int j = 0; j < m + 1; j++) {
            A[i][j] = 0;
        }
   }

    /**
    Calculate the weights using the chosen weighting function and calculate each
    term in the basis function (linear or quadratic).*/
    double (*v)[m] = malloc((np + 1) * sizeof(*v)); // basis matrix for each point 
                                                      // (plus interplated point)
    for (int n = 0; n < np + 1; n++) {
        if (n < np) {
            fill_basis(m, v[n], (coord){(p[n].x - pint.x)/h, (p[n].y - pint.y)/h});
            double r = sqrt(sq(pint.x - p[n].x) + sq(pint.y - p[n].y));
            w[n] = weighting_func(r, range);

            for (int i = 0; i < m; ++i)
                A[i][m] += w[n]*v[n][i]*u[n]; // b vector
        } else {
            fill_basis(m, v[n], (coord){0,0,0});
        }
    }

    /**
    Construct the moment matrix and vector.*/
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int n = 0; n < np; n++) {
                A[i][j] += w[n]*v[n][i]*v[n][j];
            }
        }
    }

#if 0
    fprintf(stderr, "np=%d\n", np);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m + 1; j++)
        fprintf(stderr, "A[%d][%d] = %g\n", i, j, A[i][j]);
#endif
    /**
    Solve the linear system, where A is the augmented matrix [A | b]. */
    //gauss_elim(m, m + 1, A, c);
    if (gauss_elim(m, m + 1, A, c) < 0)
        for (int n = 0; n < np; n++)
            fprintf(stderr, "|| %d p = {%g, %g} u = %g\n", n, p[n].x, p[n].y, u[n]);

    /**
    Calculate the value at the inerpolated point using the found coefficients.*/
    double uint = 0;
    for (int i = 0; i < m; i++)
        uint += v[np][i]*c[i];

    if (A) free(A); 
    if (w) free(w); 
    if (c) free(c); 
    if (v) free(v);

    return uint;
}

