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

double magnitude_coord (coord a)
{
    double mag = 0;
    foreach_dimension()
        mag += sq(a.x);
    return sqrt(mag);
}

double dot_product_angle (coord a, coord b)
{
    double product = clamp(dot_product(a, b), -1, 1);
    double amag = magnitude_coord(a);
    double bmag = magnitude_coord(b);

    return acos(product/(amag * bmag));
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
 ### Bug fix for fraction function (2D case) 

This is pulled from Tianyang's sandbox, so thank you to him.

It fixes the issue with the standard fractions() function, where it can't
reconstruct an interface if it perfectly coincides with a vertex. */

trace
void fractions_ibm (vertex scalar Phi, scalar c,
		face vector s = {0}, double val = 0.)
{
#if dimension > 1
  face vector as = automatic (s);
  
  /**
  We store the positions of the intersections of the surface with the
  edges of the cell in vector field `p`. In two dimensions, this field
  is just the transpose of the *line fractions* `s`, in 3D we need to
  allocate a new field. */
  
#if dimension == 3
  vector p[];
#else // dimension == 2
  vector p;
  p.x = as.y; p.y = as.x;
#endif
  
  /**
  ### Line fraction computation
  
  We start by computing the *line fractions* i.e. the (normalised)
  lengths of the edges of the cell within the surface. */

  foreach_edge() {

    /**
    If the values of $\Phi$ on the vertices of the edge have opposite
    signs, we know that the edge is cut by the interface. */

    if ((Phi[] - val)*(Phi[1] - val) < 0.) {

      /**
      In that case we can find an approximation of the interface position by
      simple linear interpolation. We also check the sign of one of the
      vertices to orient the interface properly. */

      p.x[] = (Phi[] - val)/(Phi[] - Phi[1]);
      if (Phi[] < val)
	p.x[] = 1. - p.x[];
    }

    /**
    If the values of $\Phi$ on the vertices of the edge have the same sign
    (or are zero), then the edge is either entirely outside or entirely
    inside the interface. We check the sign of both vertices to treat
    limit cases properly (when the interface intersects the edge exactly
    on one of the vertices). */

    else
      p.x[] = (Phi[] > val || Phi[1] > val);
  }

  /**
  ### Surface fraction computation 

  We can now compute the surface fractions. In 3D they will be
  computed for each face (in the z, x and y directions) and stored in
  the face field `s`. In 2D the surface fraction in the z-direction is
  the *volume fraction* `c`. */

#if dimension == 3

  /**
  In 3D we need to prevent boundary conditions, since this would
  impose vertex field BCs which are not (apparently) consistent for
  the edge intersection coordinates. This can probably be improved. */
  
  foreach_dimension()
    p.x.stencil.bc |= s_centered|s_face;
  
  scalar s_x = as.x, s_y = as.y, s_z = as.z;
  foreach_face(z,x,y)
#else // dimension == 2
  scalar s_z = c;
  foreach()
#endif
  {

    /**
    We first compute the normal to the interface. This can be done easily
    using the line fractions. The idea is to compute the circulation of
    the normal along the boundary $\partial\Omega$ of the fraction of the
    cell $\Omega$ inside the interface. Since this is a closed curve, we
    have
    $$
    \oint_{\partial\Omega}\mathbf{n}\;dl = 0
    $$ 
    We can further decompose the integral into its parts along the edges
    of the square and the part along the interface. For the case pictured
    above, we get for one component (and similarly for the other)
    $$
    - s_x[] + \oint_{\Phi=0}n_x\;dl = 0
    $$
    If we now define the *average normal* to the interface as
    $$
    \overline{\mathbf{n}} = \oint_{\Phi=0}\mathbf{n}\;dl
    $$
    We have in the general case
    $$
    \overline{\mathbf{n}}_x = s_x[] - s_x[1,0]
    $$
    and
    $$
    |\overline{\mathbf{n}}| = \oint_{\Phi=0}\;dl
    $$ 
    Note also that this average normal is exact in the case of a linear
    interface. */

    coord n;
    double nn = 0.;
    foreach_dimension(2) {
      n.x = p.y[] - p.y[1];
      nn += fabs(n.x);
    }
    
    /**
    If the norm is zero, the cell is full or empty and the surface fraction
    is identical to one of the line fractions. */

    if (nn == 0.)
      s_z[] = p.x[];
    else {
    
      /**
      Otherwise we are in a cell containing the interface. We first
      normalise the normal. */

      foreach_dimension(2)
	n.x /= nn;

      /**
      To find the intercept $\alpha$, we look for edges which are cut by the
      interface, find the coordinate $a$ of the intersection and use it to
      derive $\alpha$. We take the average of $\alpha$ for all intersections. */
      
      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
	foreach_dimension(2)
	  if (p.x[0,i] > 0. && p.x[0,i] < 1.) {
	    double a = sign(Phi[0,i] - val)*(p.x[0,i] - 0.5);
	    alpha += n.x*a + n.y*(i - 0.5);
	    ni++;
	  }
    else if (p.x[0,i] == 0 && p.y[i,0] == 1) { //for interface crossing grid corner
      alpha += n.x*(i - 0.5) + n.y*(i - 0.5);
      ni++;
    }
    else if (p.x[0,i] == 0 && p.y[1-i,0] == 1) {
      alpha += n.x*(0.5 - i) + n.y*(i - 0.5);
      ni++;
    }

      /**
      Once we have $\mathbf{n}$ and $\alpha$, the (linear) interface
      is fully defined and we can compute the surface fraction using
      our pre-defined function. For marginal cases, the cell is full
      or empty (*ni == 0*) and we look at the line fractions to
      decide. */

      if (ni == 0)
	s_z[] = max (p.x[], p.y[]);
      else if (ni != 4)
	s_z[] = line_area (n.x, n.y, alpha/ni);
      else {
#if dimension == 3
	s_z[] = (p.x[] + p.x[0,1] + p.y[] + p.y[1] > 2.);
#else
	s_z[] = 0.;
#endif
      }
    }
  }
  
  /**
  ### Volume fraction computation

  To compute the volume fraction in 3D, we use the same approach. */
  
#if dimension == 3
  foreach() {

    /**
    Estimation of the average normal from the surface fractions. */
       
    coord n;
    double nn = 0.;
    foreach_dimension(3) {
      n.x = as.x[] - as.x[1];
      nn += fabs(n.x);
    }
    if (nn == 0.)
      c[] = as.x[];
    else {
      foreach_dimension(3)
	n.x /= nn;

      /**
      We compute the average value of *alpha* by looking at the
      intersections of the surface with the twelve edges of the
      cube. */
      
      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
	for (int j = 0; j <= 1; j++)
	  foreach_dimension(3)
	    if (p.x[0,i,j] > 0. && p.x[0,i,j] < 1.) {
	      double a = sign(Phi[0,i,j] - val)*(p.x[0,i,j] - 0.5);
	      alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
	      ni++;
	    }

      /**
      Finally we compute the volume fraction. */

      if (ni == 0)
	c[] = as.x[];
      else if (ni < 3 || ni > 6)
	c[] = 0.; // this is important for robustness of embedded boundaries
      else
	c[] = plane_volume (n, alpha/ni);
    }
  }
#endif // dimension == 3
#else  // dimension == 1
  if (s.x.i)
    foreach_face()
      s.x[] = Phi[] > 0.;
  foreach()
    if ((Phi[] - val)*(Phi[1] - val) < 0.) {
      c[] = (Phi[] - val)/(Phi[] - Phi[1]);
      if (Phi[] < val)
	c[] = 1. - c[];
    }
    else
      c[] = (Phi[] > val || Phi[1] > val);
#endif
}

macro fraction_ibm (scalar f, double func)
{
  {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = func;
    fractions_ibm (phi, f);
  }
}

macro solid_ibm (scalar cs, face vector fs, double func)
{
  {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = func;
    fractions_ibm (phi, cs, fs);
  }
}


int facets_ibm (coord n, double alpha, coord p[2])
{
  int nx = 0;
  int ny = 0;
  coord px[2];
  coord py[2];

  for (double s = -0.5; s <= 0.5; s += 1.) {
    if (fabs (n.y) > 0.) {
      double a = (alpha - s*n.x)/n.y;
      if (a >= -0.5 && a <= 0.5) {
        px[nx].x = s;
	      px[nx++].y = a;
      }
    }
    if (fabs (n.x) > 0.) {
      double a = (alpha - s*n.y)/n.x;
      if (a >= -0.5 && a <= 0.5) {
        py[ny].y = s;
	      py[ny++].x = a;
      }
    }
  }

  if (nx == 2) {
    foreach_dimension() {
      p[0].x = px[0].x;
      p[1].x = px[1].x;
    }
    return 2;
  }
  else if (ny == 2) {
    foreach_dimension() {
      p[0].x = py[0].x;
      p[1].x = py[1].x;
    }
    return 2;
  }
  else {
    int i = 0;
    if (nx > 0) {
      assert (nx == 1);
      p[i].x = px[0].x;
      p[i].y = px[0].y;
      i++;
    }
    if (ny > 0) {
      assert (ny == 1);
      p[i].x = py[0].x;
      p[i].y = py[0].y;
      i++;
    }
    assert (i <= 2);
    return i;
  }
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
This like the normal statsf function, but without dv(). It also accepts an
input parameter which tells to check for cells where f != val. For example,
f != val = 0. val = nodata by default.

The *statsf()* function returns the minimum, maximum, volume sum,
standard deviation and volume for field *f*. */

typedef struct {
  double min, max, sum, stddev, volume;
} stats2;

stats2 statsf2 (scalar f, double val = nodata)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	  reduction(max:max) reduction(min:min)) 
    if (f[] != nodata) {
      volume += 1;
      sum    += f[];
      sum2   += sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats2 s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
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
    *rank = -1;

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
      else if (allocated(0) && is_boundary(cell)) {
        if (rank)
            *rank = -2;
        //fprintf(stderr, "is boundary cell! (%g, %g, %g)\n", xp, yp, zp);
      }
    }
    else
      break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

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



