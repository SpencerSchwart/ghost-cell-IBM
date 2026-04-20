/**
Header file for enforcing a pinned contact line using contact.h */

#include "contact.h"

scalar f0[]; // f at previous time step

vector h[];
double theta0 = 0.5*pi; // default ca @t=0 is 90 degrees.

int contact_boundary = left;

double dot_product (coord a, coord b)
{
    double product = 0;
    foreach_dimension()
        product += a.x*b.x;
    return product;
}

double cross_product_2d (coord a, coord b)
{
    return (a.x*b.y) - (a.y*b.x);
}


void normalize_sum (coord * n)
{
    double norm = 0;
    foreach_dimension()
        norm += fabs(n->x);
    foreach_dimension()
        n->x /= norm + SEPS;
}

int approx_equal_double (double a, double b, double TOL = 1e-10)
{
    return fabs(a - b) <= TOL;
}

/**
generic root solver using Brent's method */
int rsolver_brent (double* result, double a, double b, const void* data, double (*func)(const void*, double),
                   double tolerance = 1e-10, int maxitr = 50, int warning = 1)
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

/**
line_intersect returns the x (resp. y) coordinate along a line given y (resp. x).
The return value is in the cell's local coordinate system, i.e. [-0.5,-0.5]x[0.5,0.5]. */
double line_intersect (double alpha, coord n, double x = HUGE, double y = HUGE)
{
    double inter = 0;
    if (x == HUGE && y != HUGE)
        inter = alpha/n.x - (n.y/(n.x+SEPS))*y;
    else
        inter = alpha/n.y - (n.x/(n.y+SEPS))*x;
    return inter;
}

/**
calculate the desired normal based on the contact angle */
static inline coord normal_contact (coord ns, coord nf, double angle)
{
    coord n;
    if (- ns.x * nf.y + ns.y * nf.x >= 0) { // 2D cross product
        n.x = - ns.x * cos(angle) + ns.y * sin(angle);
        n.y = - ns.x * sin(angle) - ns.y * cos(angle);
    }
    else {
        n.x = - ns.x * cos(angle) - ns.y * sin(angle);
        n.y =   ns.x * sin(angle) - ns.y * cos(angle);
    }
    return n;
}

/**
for root solver */
typedef struct hysteresis_error
{
    double c, alphas;
    coord pint, ns;
    int csign; // for orientating the interface correctly
} hysteresis_error;

double get_hysteresis_error (const void* data, double theta)
{
    const hysteresis_error tcell = (*(hysteresis_error*)data);
    coord ns = tcell.ns, nf;

    if (tcell.csign > 0) {
        nf.x = ns.x*cos(theta) - ns.y*sin(theta);
        nf.y = ns.x*sin(theta) + ns.y*cos(theta);
    }
    else {
        nf.x =  ns.x*cos(theta) + ns.y*sin(theta);
        nf.y = -ns.x*sin(theta) + ns.y*cos(theta);       
    }
    double alphaf = tcell.alphas + dot_product(nf, tcell.pint) - dot_product(ns, tcell.pint);

    double fa = rectangle_fraction (nf, alphaf, (coord){-0.5,-0.5}, (coord){0.5,0.5});
    return fa - tcell.c;
}

/**
returns the angle (in radians) of the contact line that keeps the interface pinned at point = pint
while c updates from c^n -> c^n+1. */
double get_hysteresis_angle (double c, coord pint, coord ns, double alphas, int csign)
{
    hysteresis_error tcell = {c, alphas, pint, ns, csign};

    double theta = 0;

    double theta_min = 0, theta_max = M_PI; // range

    rsolver_brent (&theta, theta_min, theta_max, &tcell, get_hysteresis_error, tolerance = 1e-10);
    return csign? M_PI - theta: theta;
}

/**
returns outward pointing normal vector of a boundary. */
coord get_boundary_normal(int boundary)
{
  if (boundary == left)
    return (coord){-1,0};
  if (boundary == right)
    return (coord){1,0};
  if (boundary == top)
    return (coord){0,1};
  if (boundary == bottom)
    return (coord){0,-1};
  return (coord){0,0};
}

/**
returns the intercept (alpha) of the boundary */
double get_boundary_intercept(int boundary)
{
  if (boundary == left || boundary == bottom)
    return 0.5;
  if (boundary == right || boundary == top)
    return -0.5;
  return 0;
}

void update_angle()
{
  /** boundary variables */
  coord ns = get_boundary_normal(contact_boundary);
  double alphas = get_boundary_intercept(contact_boundary);

  /** find contact line cell */
  foreach_boundary(contact_boundary) {
    if (f[] > 0 && f[] < 1 && f0[] > 0 && f0[] < 1) { // potential contact line cell

      coord nf0 = interface_normal(point, f0);
      coord nc0 = normal_contact (ns, nf0, theta0); // desired normal
      normalize_sum(&nc0);

      double alphaf0 = plane_alpha(f0[], nc0); // desired intercept

      /** calculate contact point coordinates */
      coord pcl0;
      if (contact_boundary == bottom || contact_boundary == top)
        pcl0 = (coord){line_intersect(alphaf0, nc0, y = sign2(ns.y)*0.5), y = sign2(ns.y)*0.5};
      else
        pcl0 = (coord){sign2(ns.x)*0.5, line_intersect(alphaf0, nc0, x = sign2(ns.x)*0.5)};

      /** contact line point must lie within cell boundary */
      if (fabs(pcl0.x) <= 0.5 && fabs(pcl0.y) <= 0.5) {
        int csign = sign2(cross_product_2d(ns, nc0)); // orient the interface correctly
        theta0 = get_hysteresis_angle(f[], pcl0, ns, alphas, csign);
        if (theta0 <= 0)
            theta0 = 0.1;
      }
    }
  }

  /** update contact angle boundary condition */
  h.t[contact_boundary] = contact_angle (theta0);
}

event defaults (i = 0)
{
  f0.height = h;
  f.height = h;

  h.t[contact_boundary] = contact_angle (theta0);
}

/**
update the contact angle (after f is updated in the vof event) */
event tracer_advection(i++)
{
  if (i == 0) {
    foreach()
      f0[] = f[];
  }

  update_angle();
}

/**
update f0 */
event end_timestep (i++)
{
  foreach()
    f0[] = f[];
}
