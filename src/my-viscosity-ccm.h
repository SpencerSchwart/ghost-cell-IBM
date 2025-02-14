#include "my-poisson.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
};

#define lambda ((coord){0.,0.,0.})




#undef face_gradient_x
#define face_gradient_x(a,i)					\
  (a.third && fs.x[i] < 1. && fs.x[i] > 0. ?			\
   embed_face_gradient_x (point, a, i) :			\
   (a[i] - a[i-1])/Delta)

#undef face_gradient_y
#define face_gradient_y(a,i)					\
  (a.third && fs.y[0,i] < 1. && fs.y[0,i] > 0. ?		\
   embed_face_gradient_y (point, a, i) :			\
   (a[0,i] - a[0,i-1])/Delta)



static void relax_diffusion (scalar * a, scalar * b, int l, void * data)
{
    struct Viscosity * p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]);

    foreach_level_or_leaf (l) {
        double avgmu = 0.;
        foreach_dimension()
              avgmu += ibmf.x[]*mu.x[] + ibmf.x[1]*mu.x[1];
        avgmu = dt * avgmu + SEPS;
        
        foreach_dimension() {
            double c = 0.;
            double d = ibm_flux_x (point, u.x, mu, &c);
            scalar s = u.x;
            double a = 0.;

            foreach_dimension()
                a += ibmf.x[1]*mu.x[1]*s[1] + ibmf.x[]*mu.x[]*s[-1];

            u.x[] = (dt * a + (r.x[] - dt*c) * sq(Delta)) /
                    (sq(Delta) * (rho[] + lambda.x + dt*d) + avgmu); 
        }
    }
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
#if TREE
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
        }
    }
#else
    foreach (reduction(max:maxres), nowarning) {
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

    restriction ({mu, rho, cm, ibmf});
    struct Viscosity p = { mu, rho, dt };
    return mg_solve ((scalar *){u}, (scalar *){r},
                     residual_diffusion, relax_diffusion, &p, 
                     nrelax, res, minlevel = 1, 
                     tolerance = TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE);
}

