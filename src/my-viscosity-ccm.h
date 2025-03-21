#include "my-poisson.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
  double (* ibm_flux_x) (Point, scalar, vector, double *);
  double (* ibm_flux_y) (Point, scalar, vector, double *);
};

#define lambda ((coord){0.,0.,0.})

#if 0
#define face_condition(ibmf, ibm)						\
  (ibmf.x[i,j] > 0.5 && ibmf.y[i,j + (j < 0)] && ibmf.y[i-1,j + (j < 0)] &&	\
   ibm[i,j] && ibm[i-1,j])

foreach_dimension()
static inline double ibm_face_gradient_x (Point point, scalar a, int i)
{
  int j = sign(ibmf.x[i,1] - ibmf.x[i,-1]);
  assert (ibm[i] && ibm[i-1]);
  if (face_condition (ibmf, ibm))
    return ((1. + ibmf.x[i])*(a[i] - a[i-1]) +
	    (1. - ibmf.x[i])*(a[i,j] - a[i-1,j]))/(2.*Delta);
  return (a[i] - a[i-1])/Delta;
}

#undef face_gradient_x
#define face_gradient_x(a,i)					\
  (ibmf.x[i] < 1. && ibmf.x[i] > 0. ?			\
   ibm_face_gradient_x (point, a, i) :			\
   (a[i] - a[i-1])/Delta)

#undef face_gradient_y
#define face_gradient_y(a,i)					\
  (ibmf.y[0,i] < 1. && ibmf.y[0,i] > 0. ?		\
   ibm_face_gradient_y (point, a, i) :			\
   (a[0,i] - a[0,i-1])/Delta)
#endif


static void relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
    struct Viscosity * p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]);

    double (* ibm_flux_x) (Point, scalar, vector, double *) = p->ibm_flux_x;
    double (* ibm_flux_y) (Point, scalar, vector, double *) = p->ibm_flux_y;
    foreach_level_or_leaf (l, nowarning) {
        double avgmu = 0.;
        foreach_dimension()
            avgmu += ibmf.x[]*mu.x[] + ibmf.x[1]*mu.x[1];
        avgmu = dt*avgmu + SEPS;
        foreach_dimension() {
            double c = 0.;
            double d = ibm_flux_x ? ibm_flux_x (point, u.x, mu, &c) : 0.;
            scalar s = u.x;
            double a = 0.;
            foreach_dimension()
            	a += ibmf.x[1]*mu.x[1]*s[1] + ibmf.x[]*mu.x[]*s[-1];
            u.x[] = (dt*a + (r.x[] - dt*c)*sq(Delta))/
                    (sq(Delta)*(rho[] + lambda.x + dt*d) + avgmu);
        }
    }
  

#if 0
    foreach_level_or_leaf (l, nowarning) {
        double avgmu = 0.;
        foreach_dimension()
              avgmu += ibmf.x[]*mu.x[] + ibmf.x[1]*mu.x[1];
        avgmu = dt * avgmu + SEPS;
        
        coord ind = {0,1};
        foreach_dimension() {
            double c = 0.;
            double d = ibm_flux_x (point, u.x, mu, &c);
            scalar s = u.x;
            double a = 0.;

            foreach_dimension()
                a += ibmf.x[1]*mu.x[1]*s[1] + ibmf.x[]*mu.x[]*s[-1];

            u.x[] = (dt * a + (r.x[] - dt*c) * sq(Delta)) /
                    (sq(Delta) * (rho[] + lambda.x + dt*d) + avgmu); 

            if (fabs(a) > LIMIT) {
                fprintf(stderr, "\n| WARNING in viscosity relax: u_%g = %g in (%g, %g) exceeds %g\n",
                                  ind.x, u.x[], x, y, LIMIT);
                fprintf(stderr, "| lvl=%d d=%g c=%g a=%g r_%g=%g rho=%g avgmu=%g ibm=%g\n", 
                                  l, d, c, a, ind.x, r.x[], rho[], avgmu, ibm[]);
            }
        }
    }
#endif    
}

static double residual_diffusion (scalar * a, scalar * b, scalar * resl,
                                  void * data)
{
    struct Viscosity * p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
    double maxres = 0.;
    foreach_dimension() {
      scalar s = u.x;
        face vector g[];
        foreach_face()
            g.x[] = ibmf.x[]*mu.x[]*face_gradient_x (s, 0);
        foreach (reduction(max:maxres), nowarning) {
            double a = 0.;
            foreach_dimension()
            	a += g.x[] - g.x[1];
            res.x[] = r.x[] - (rho[] + lambda.x)*u.x[] - dt*a/Delta;
            if (ibm_flux_x) {
            	double c, d = ibm_flux_x (point, u.x, mu, &c);
            	res.x[] -= dt*(c + d*u.x[]);
              }
            if (fabs (res.x[]) > maxres)
            	maxres = fabs (res.x[]);
        }
    }
#if 0
#if TREE
    coord ind = {0,1};
    foreach_dimension() {
        scalar s = u.x;
        face vector g[];

        foreach_face() {
            g.x[] = ibmf.x[] * mu.x[] * face_gradient_x (s, 0);
        }

        foreach (reduction(max:maxres), nowarning) {
            double a = 0.;
            foreach_dimension()
                a += g.x[] - g.x[1];
            res.x[] = r.x[] - (rho[] + lambda.x) * u.x[] - dt * a / Delta;

            double c, d = ibm_flux_x (point, u.x, mu, &c);
            res.x[] -= dt * (c + d*u.x[]);

            if (fabs (res.x[]) > maxres)
                maxres = fabs (res.x[]);

            if (fabs(res.x[]) > LIMIT) {
                fprintf(stderr, "\n| WARNING in viscosity residual: res_%g = %g in (%g, %g) exceeds %g\n",
                                  ind.x, res.x[], x, y, LIMIT);
                fprintf(stderr, "| d=%g c=%g a=%g r_%g=%g rho=%g g[]=%g ibm=%g\n", 
                                  d, c, a, ind.x, r.x[], rho[], g.x[], ibm[]);
            }
        }
    }
#else
    foreach (reduction(max:maxres), nowarning)
        foreach_dimension() {
            scalar s = u.x;
            double a = 0.;
            foreach_dimension()
                a += mu.x[0]*face_gradient_x (s, 0) - mu.x[1]*face_gradient_x (s, 1);
            res.x[] = r.x[] - (rho[]/(ibm[] + SEPS) + lambda.x) * u.x[] - dt * a / Delta;

            double c, d = ibm_flux_x (point, u.x, mu, &c);
            res.x[] -= dt * (c + d*u.x[]);

            if (fabs (res.x[]) > maxres)
                maxres = fabs (res.x[]);
        }
#endif
#endif
    return maxres;
}

#undef lambda

double TOLERANCE_MU = 0.;

trace
mgstats viscosity (vector u, face vector mu, scalar rho, double dt,
                   int nrelax = 4, scalar * res = NULL)
{
    vector r[];
    foreach() {
        foreach_dimension() {
            r.x[] = rho[]/(ibm[] + SEPS) * u.x[];
        }
    }

    restriction ({mu, rho, cm, ibmf, ibm});
    struct Viscosity p = { mu, rho, dt };
    #if 1
    p.ibm_flux_x = ibm_flux_x;
    p.ibm_flux_y = ibm_flux_y;
    #else
    p.ibm_flux_x = NULL;
    p.ibm_flux_y = NULL;
    #endif
    return mg_solve ((scalar *){u}, (scalar *){r},
                     residual_diffusion, relax_diffusion, &p, 
                     nrelax, res, minlevel = 1, 
                     tolerance = TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE);
}

#if 0
#undef face_gradient_x
#define face_gradient_x(a,i) ((a[i] - a[i-1])/Delta)
#undef face_gradient_y
#define face_gradient_y(a,i) ((a[0,i] - a[0,i-1])/Delta)
#undef face_gradient_z
#define face_gradient_z(a,i) ((a[0,0,i] - a[0,0,i-1])/Delta)
#endif
