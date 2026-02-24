/**
Header file containing helper functions for IBM. */

#define VTOL 1e-11
#define INT_TOL 1e-7    // tolerance used for volume fraction fields (interface tolerance)

#define distance(a,b) sqrt(sq(a) + sq(b))
#define distance3D(a,b,c) sqrt(sq(a) + sq(b) + sq(c))

#define on_interface(a) (a[] > 0+INT_TOL && a[] < 1-INT_TOL)
//#define on_interface(a,_TOL) (a[] > 0+_TOL && a[] < 1.-_TOL)

#define is_mostly_solid(a, i) (a[i] > 0+INT_TOL && a[i] <= 0.5)
#define is_fresh_cell(a0, a) (a0[] <= 0.5 && a[] > 0.5)

double cross_product_2d (coord a, coord b)
{
    return (a.x*b.y) - (a.y*b.x);
}

coord cross_product (coord a, coord b)
{
    coord c = {
        (a.y*b.z) - (a.z*b.y) ,
      -((a.x*b.z) - (b.x*a.z)),
        (a.x*b.y) - (a.y*b.x)
    };
    return c;
}

double determinant (coord a, coord b)
{
    coord c = cross_product(a,b);
    return distance3D(c.x,c.y,c.z); // area
}


int approx_equal (coord p1, coord p2, double TOL = VTOL)
{
    return fabs(p1.x - p2.x) <= TOL && fabs(p1.y - p2.y) <= TOL && fabs(p1.z - p2.z) <= TOL;
}

int approx_equal_double (double a, double b, double TOL = VTOL)
{
    return fabs(a - b) <= TOL;
}


// normalize() but with SEPS in the denominator
void normalize2 (coord * n)
{
    double norm = 0;
    foreach_dimension()
        norm += sq(n->x);
    norm = sqrt(norm);
    foreach_dimension()
        n->x /= norm + SEPS;
}

void normalize_sum (coord * n)
{
    double norm = 0;
    foreach_dimension()
        norm += fabs(n->x);
    foreach_dimension()
        n->x /= norm + SEPS;
}

double dot_product_norm (coord a, coord b)
{
    normalize2(&a); normalize2(&b);

    double product = 0;
    foreach_dimension()
        product += a.x*b.x;

    return product;
}

double dot_product (coord a, coord b)
{
    double product = 0;
    foreach_dimension()
        product += a.x*b.x;

    return product;
}

double dot_product_angle (coord a, coord b)
{
    double product = clamp(dot_product(a, b), -1, 1);
    assert(fabs(product) <= 1);
    return acos(dot_product(a, b));
}

/*
normal_and_tangents accepts a normal vector and fill t1 and t2 with the 
corresponding tangent vector(s) (one in 2D and two in 3D), all normalized.
*/
void normal_and_tangents (coord * n, coord * t1, coord * t2)
{
    normalize2(n);
#if dimension == 2

    coord t1_tmp = {-n->y, n->x};
    *t1 = t1_tmp;
    *t2 = (coord){0,0,0};

#else // dimension == 3

    coord a = {0,0,1};
    if ((!fabs(dot_product(*n, a))) < 0.9) {
        a = (coord){1,0,0};
    }
    coord t1_tmp = cross_product(*n, a);
    double det = distance3D(t1_tmp.x, t1_tmp.y, t1_tmp.z); // determinant

    assert(fabs(det) > 1e-15);

    foreach_dimension()
        t1->x = t1_tmp.x/det;
    
    *t2 = cross_product(*n, *t1);

#endif // dimension == 3
}

static inline double vertex_average (Point point, scalar s)
{
#if dimension == 2
    return (4.*s[] + 
	        2.*(s[0,1] + s[0,-1] + s[1,0] + s[-1,0]) +
    	    s[-1,-1] + s[1,-1] + s[1,1] + s[-1,1])/16.;
#else
    return (8.*s[] +
	        4.*(s[-1] + s[1] + s[0,1] + s[0,-1] + s[0,0,1] + s[0,0,-1]) +
	        2.*(s[-1,1] + s[-1,0,1] + s[-1,0,-1] + s[-1,-1] + 
    		s[0,1,1] + s[0,1,-1] + s[0,-1,1] + s[0,-1,-1] +
	    	s[1,1] + s[1,0,1] + s[1,-1] + s[1,0,-1]) +
	        s[1,-1,1] + s[-1,1,1] + s[-1,1,-1] + s[1,1,1] +
    	    s[1,1,-1] + s[-1,-1,-1] + s[1,-1,-1] + s[-1,-1,1])/64.;
#endif
}


/*
local_to_global fills ax, ay, and az with the global coordinates of a point, p,
given in a local coordinate system.
*/

int local_to_global (Point point, coord p, double* ax, double* ay, double* az)
{
    *ax = x + p.x * Delta;
    *ay=  y + p.y * Delta;
    *az = z + p.z * Delta;

    return 0;
}

/**
Same as above, but returns the values in a coord type instead of individual double variables.*/
coord local_to_global_coord (Point point, coord p1)
{
    coord c = {x,y,z}, p2;
    foreach_dimension()
        p2.x = c.x + p1.x*Delta;
    return p2;
}



/*
copy_coord is used to fill three variables with the coresponding components of p.
This is used to avoid having the .x or _x indicies being automatically changed
within a foreach_dimension()
*/

int copy_coord (coord p, double* ax, double* ay, double* az)
{
    *ax = p.x;
    *ay = p.y;
    *az = p.z;

    return 1;
}

/*
gauss_elim performs *in place* transformation to the provided augmented matrix 
(meaning it is changed w/o making a copy) and fills coeff with the solved linear system.

***Courtesy of ChatGPT***

TODO: Extend to handle higher-order interpolation schemes, i.e. larger matrices.
*/

//void gauss_elim(int m, int n, double matrix[m][n], double sol[m])
int gauss_elim(int m, int n, double matrix[m][n], double sol[m])
{
    // Forward elimination
    for (int i = 0; i < m; i++) {

        // 1. Partial pivot: find row with largest pivot in column i
        int max_row = i;
        for (int r = i + 1; r < m; r++) {
            if (fabs(matrix[r][i]) > fabs(matrix[max_row][i])) {
                max_row = r;
            }
        }

        // 2. Swap current row i with max_row if needed
        if (max_row != i) {
            for (int c = 0; c < n; c++) {
                double temp = matrix[i][c];
                matrix[i][c]  = matrix[max_row][c];
                matrix[max_row][c] = temp;
            }
        }

        // 3. Make sure our pivot is non‐zero (or not too close to zero)
        if (fabs(matrix[i][i]) < 1e-14) {
            fprintf(stderr, "ERROR: Pivot is zero (matrix is singular or nearly singular)\n");
            return -1;
        }

        // 4. Eliminate all rows below row i
        for (int r = i + 1; r < m; r++) {
            double factor = matrix[r][i] / matrix[i][i];

            for (int c = i; c < n; c++) {
                matrix[r][c] -= factor * matrix[i][c];
            }
        }
    }

    // Back‐substitution
    for (int i = m - 1; i >= 0; i--) {
        // Start with the RHS of the augmented matrix
        sol[i] = matrix[i][n - 1];

        // Subtract the known terms from columns to the right
        for (int c = i + 1; c < m; c++) {
            sol[i] -= matrix[i][c] * sol[c];
        }

        // Divide by the diagonal element
        sol[i] /= matrix[i][i];
    }
    return 1;
}

/**
reconstruction_ibm accepts an additional face vector field that represents the
face solid volume fraction (fs) to be used when calculating the normal, n.

TODO: what if interface perfectly cuts cell face? use interfacial() instead? must
      be as cheap as possible.
*/

void reconstruction_cs (const scalar c, const face vector cf, vector n, scalar alpha)
{
    foreach() {
        if (c[] <= 0. || c[] >= 1.) {
            alpha[] = 0.;
            foreach_dimension()
                n.x[] = 0.;
        }
        else {
            coord m = facet_normal (point, c, cf);
            foreach_dimension()
                n.x[] = m.x;
            alpha[] = plane_alpha(c[], m);
        }
    }

#if TREE
    foreach_dimension()
        n.x.refine = n.x.prolongation = refine_injection;

    alpha.n = n;
    alpha.refine = alpha.prolongation = alpha_refine;
#endif
}

/**
generic root solver using the bisection method */
int rsolver_bisection (double* a, double amin, double amax, const void* data, double (*func)(const void*, double),
                       double tolerance = BI_TOL, int maxitr = 50)
{
    double errmax = (*func)(data, amax);
    if (fabs(errmax) < tolerance) {
        *a = amax;
        return 0;
    }

    double b = 0, error = HUGE;
    int itr = 0;
    while (fabs(error) > tolerance && itr < maxitr) {
        b = (amin + amax)/2.;   // bisection
        error = (*func)(data, b);
        if (sign2(error) == sign2(errmax)) {
            amax = b;
            errmax = error;
        }
        else
            amin = b;
        itr++;
    }
        
    if (itr == maxitr) {
        fprintf(stderr, "WARNING: alpha  solver does not converge after"
                        " maximum iteration (%d), error = %g\n", maxitr, error);
    }
    *a = b;
    return itr;
}

/**
generic root solver using Brent's method */
int rsolver_brent (double* result, double a, double b, const void* data, double (*func)(const void*, double),
                   double tolerance = BI_TOL, int maxitr = 50, int warning = 1)
{
    double fa = (*func)(data, a);
    if (fabs(fa) < tolerance) {
        *result = a;
        return 0;
    }

    double fb = (*func)(data, b);
    if (fabs(fb) < tolerance) {
        *result = b;
        return 0;
    }

    if (fa * fb >= 0) {
        if (warning)
            fprintf(stderr, "WARNING: range in Brent's solver does not contain root!\n");
        return -1;
    }

    if (fabs(fa) < fabs(fb)) {
        swap(double, a, b);
        swap(double, fa, fb);
    }

    double c = a, fc = fa, d = b - a; 
    int iter = 0, mflag = 1;
    while (fabs(b - a) > tolerance && fabs(fb) > tolerance && fabs(fa) > tolerance && iter < maxitr) {

        double s;
        if (!approx_equal_double(fa, fc) && !approx_equal_double(fb, fc)) {
            s = (a*fb*fc)/((fa-fb)*(fa-fc)) + 
                (b*fa*fc)/((fb-fa)*(fb-fc)) + 
                (c*fa*fb)/((fc-fa)*(fc-fb));
        }
        else // secant
            s = b - fb*(b-a)/(fb-fa + SEPS);

        int cond1 = s <= fmin(a, b) || s >= fmax(a, b);
        //int cond1 = (s < (3*a + b)/4.) || (s > b); // assuming a < b and |fa| > |fb|
        int cond2 =  mflag && fabs(s - b) >= fabs(b-c)*0.5;
        int cond3 = !mflag && fabs(s - b) >= fabs(c-d)*0.5;
        int cond4 =  mflag && fabs(b - c) < tolerance;
        int cond5 = !mflag && fabs(c - d) < tolerance;
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = (a + b) * 0.5;
            mflag = 1;
        }
        else
            mflag = 0;

        double fs = (*func)(data, s);
        d = c;        
        c = b; fc = fb;

        if (fa*fs < 0)
            b = s, fb = fs;
        else
            a = s, fa = fs;

        if (fabs(fa) < fabs(fb)) {
            swap(double, a, b);
            swap(double, fa, fb);
        }
        if (fabs(fb) < tolerance || fabs(fa) < tolerance)
            break;

        ++iter;
    }

    if (iter == maxitr && warning) 
        fprintf(stderr, "WARNING: alpha  solver does not converge after"
                        " maximum iteration (%d), error = %g\n", maxitr, fb);

    *result = (fabs(fa) < fabs(fb)) ? a : b;
    return iter;
}


