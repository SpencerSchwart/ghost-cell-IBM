/**
Header file containing helper functions for IBM. */

#define VTOL 1e-11
#define INT_TOL 1e-6    // tolerance used for volume fraction fields (interface tolerance)

#define distance(a,b) sqrt(sq(a) + sq(b))
#define distance3D(a,b,c) sqrt(sq(a) + sq(b) + sq(c))

#if dimension == 2
#define distance_coord(a,b) (sqrt(sq(a.x - b.x) + sq(a.y - b.y)))
#else
#define distance_coord(a,b) (sqrt(sq(a.x - b.x) + sq(a.y - b.y) + sq(a.z - b.z)))
#endif

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

/**
Calculate the determinant of a 3x3 matrix, of which the row vectors a1, 
a2, and a3 are given as input. */
double determinant_rows (coord a1, coord a2, coord a3)
{
    return a1.x*a2.y*a3.z + a1.y*a2.z*a3.x + a1.z*a2.x*a3.y - 
           a1.x*a2.z*a3.y - a1.y*a2.x*a3.z - a1.z*a2.y*a3.x;
}

/**
Calculate the determinant of a 3x3 matrix, of which the columns vectors a1, 
a2, and a3 are given as input. */
static inline double determinant_cols (coord a1, coord a2, coord a3)
{
    return a1.x*a2.y*a3.z + a2.x*a3.y*a1.z + a3.x*a1.y*a2.z - 
           a1.x*a3.y*a2.z - a2.x*a1.y*a3.z - a3.x*a2.y*a1.z;
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

double magnitude_coord (coord a)
{
    double mag = 0;
    foreach_dimension()
        mag += sq(a.x);
    return sqrt(mag);
}

coord subtract_coord (coord a, coord b)
{
    coord c;
    foreach_dimension()
        c.x = a.x - b.x;
    return c; 
}

bool equals_coord (coord a, coord b)
{
#if dimension == 2
    return a.x == b.x && a.y == b.y;
#else
    return a.x == b.x && a.y == b.y && a.z == b.z;
#endif
}

bool empty_coord (coord a)
{
#if dimension == 2
    return !a.x && !a.y;
#else
    return !a.x && !a.y && !a.z;
#endif
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

void vorticity2 (const vector u, scalar omega)
{
  foreach() {
    if (cm[]) {
      omega[] = ((fm.x[1] - fm.x[])*u.y[] +
             fm.x[1]*u.y[1] - fm.x[]*u.y[-1] -
             (fm.y[0,1] - fm.y[])*u.x[] +
             fm.y[]*u.x[0,-1] - fm.y[0,1]*u.x[0,1])/(2.*(cm[] + SEPS)*Delta);
    }
  }
}

/**
## Simple field statistics 

The *normf()* function returns the (volume) average, RMS norm, max
norm and volume for field *f*. */

typedef struct {
  double avg, rms, l2, max, volume;
} norm2;

norm2 normf2 (scalar f)
{
  double avg = 0., rms = 0., l2 = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
	  reduction(+:rms) reduction(+:l2) reduction(+:volume)) 
    if (f[] != nodata && dv() > 0.) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += dv();
      avg    += dv()*v;
      rms    += dv()*sq(v);
      l2     += sq(v);
    }
  norm2 n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.l2 = l2;
  n.max = max;
  n.volume = volume;
  return n;
}

/**
The foreach_direct_neighbor macro iterates over the cells that share a face
with the given cell in the stencil. */
#if dimension == 2
macro foreach_direct_neighbor (int self = 0, Point point = point, break = (_k = _l = 2)) {
  {
    const int _ig = point.i, _jg = point.j;
    int ig = 0, jg = 0;

    for (int _k = -1; _k <= 1; _k++) {
      point.i = _ig + _k;
      for (int _l = -1; _l <= 1; _l++) {
	    point.j = _jg + _l;

        bool allow = !self? !(_l == 0 && _k == 0) && (abs(_l) != abs(_k)):
                             (_l == 0 && _k == 0) || (abs(_l) != abs(_k));
        if (allow) {
	      POINT_VARIABLES();
	      {...}
        }
      }
    }
    point.i = _ig; point.j = _jg;
  }
}
#else
macro foreach_direct_neighbor (int self = 0, Point point = point, break = (_l = _m = _n = 2)) {
  {
    const int _i = point.i, _j = point.j, _k = point.k;
    int ig = 0, jg = 0, kg = 0;

    for (int _l = -1; _l <= 1; _l++) {
      point.i = _i + _l;
      for (int _m = -1; _m <= 1; _m++) {
        point.j = _j + _m;
        for (int _n = -1; _n <= 1; _n++) {
          point.k = _k + _n;

          bool check = false;
          if (abs(_l) > 0)
            check = _m == 0 && _n == 0;
          else if (abs(_m) > 0)
            check = _l == 0 && _n == 0;
          else if (abs(_n) > 0)
            check = _l == 0 && _m == 0;
            
          bool allow = !self? !(_l == 0 && _m == 0 && _n == 0) && check:
                               (_l == 0 && _m == 0 && _n == 0) || check;
          if (allow) {
            POINT_VARIABLES();
	        {...}
          }
        }
      }
    }
    point.i = _i; point.j = _j; point.k = _k;
  }
}
#endif

#if dimension == 2
macro foreach_near_neighbor (int self = 0, Point point = point, break = (_k = _l = _nn + 1)) {
  {
    const int _nn = 1;
    const int _ig = point.i, _jg = point.j;
    int ig = 0, jg = 0;
    for (int _k = - _nn; _k <= _nn; _k++) {
      point.i = _ig + _k;
      for (int _l = - _nn; _l <= _nn; _l++) {
	    point.j = _jg + _l;
        bool allow_center = !self? !(_l == 0 && _k == 0): true;
        if (allow_center) 
        {
	      POINT_VARIABLES();
	      {...}
        }
      }
    }
    point.i = _ig; point.j = _jg;
  }
}
#else
macro foreach_near_neighbor (int self = 0, Point point = point, break = (_l = _m = _n = 2)) {
  {
    const int _i = point.i, _j = point.j, _k = point.k;
    int ig = 0, jg = 0, kg = 0;

    for (int _l = -1; _l <= 1; _l++) {
      point.i = _i + _l;
      for (int _m = -1; _m <= 1; _m++) {
        point.j = _j + _m;
        for (int _n = -1; _n <= 1; _n++) {
          point.k = _k + _n;
            
          bool allow = !self? !(_l == 0 && _m == 0 && _n == 0): true;
          if (allow) {
            POINT_VARIABLES();
	        {...}
          }
        }
      }
    }
    point.i = _i; point.j = _j; point.k = _k;
  }
}
#endif

#if dimension == 2
macro foreach_diagonal_neighbor (int self = 0, Point point = point, break = (_k = _l = _nn + 1)) {
  {
    const int _nn = 1;
    const int _ig = point.i, _jg = point.j;
    int ig = 0, jg = 0;
    for (int _k = - _nn; _k <= _nn; _k++) {
      point.i = _ig + _k;
      for (int _l = - _nn; _l <= _nn; _l++) {
	    point.j = _jg + _l;
        bool allow_center = !self? !(_l == 0 && _k == 0)  && abs(_l) == abs(_k): (abs(_l) == abs(_k)) || (_l == 0 && _k == 0);
        if (allow_center) 
        {
	      POINT_VARIABLES();
	      {...}
        }
      }
    }
    point.i = _ig; point.j = _jg;
  }
}
#else

#endif

Point locate_ibm (double xp = 0., double yp = 0., double zp = 0., int * rank = 0)
{
  if (rank)
    *rank = -2;

  for (int l = depth(); l >= 0; l--) {
    Point point = {0};
    point.level = l;
    int n = 1 << point.level;
    point.i = (xp - X0)/L0*n + GHOSTS;
#if dimension >= 2
    point.j = (yp - Y0)/L0*n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (zp - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS
#if dimension >= 2
        && point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
        && point.k >= 0 && point.k < n + 2*GHOSTS
#endif
	) {
      if (allocated(0) && is_local(cell) && is_leaf(cell)) {
        if (rank) 
            *rank = cell.pid;
    	return point;
      }
      else if (allocated(0) && !is_local(cell) && is_leaf(cell)) {
        if (rank)
            *rank = cell.pid;
        point.level = -1;
        return point;
      }
    }
    else
      break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

#if 0
macro2 foreach_image_point_stencil (double xp, double yp, double zp, char flags, Reduce reductions)
{
  foreach_stencil (flags, reductions)
    {...}
}
#endif

macro2 foreach_image_point (double _x = 0., double _y = 0., double _z = 0., int * rank = NULL,
            		        char flags = 0, Reduce reductions = None)
{
  {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    coord _p = { _x, _y, _z };
    Point point = locate_ibm (_p.x, _p.y, _p.z, rank); // fixme
    if (point.level >= 0)
      {...}
  }
}

