#include "my-poisson.h"

struct Viscosity {
  face vector mu;
  scalar rho;
  double dt;
};

#define lambda ((coord){0.,0.,0.})

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
              avgmu += mu.x[] + mu.x[1];
        avgmu = dt * avgmu + SEPS;

        foreach_dimension() {
            scalar s = u.x;
            double a = 0.;
            foreach_dimension()
                a += mu.x[1] * s[1] + mu.x[] * s[-1];
            u.x[] = cm[]*(dt * a + r.x[] * sq(Delta)) /
                         (sq(Delta) * (rho[]/(ibm[] + SEPS) + lambda.x) + avgmu); 
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
    coord ind = {0,1,2};
    foreach_dimension() {
        scalar s = u.x;
        face vector g[];

        foreach_face() {
            g.x[] = mu.x[] * face_gradient_x (s, 0);
        }

        foreach (reduction(max:maxres), nowarning) {
            double a = 0.;
            foreach_dimension()
                a += g.x[] - g.x[1];
            res.x[] = r.x[] - (rho[]/(ibm[] + SEPS) + lambda.x) * u.x[] - dt * a / Delta;
            if (ibm[] <= 0.5)
                res.x[] = 0;

            if (fabs (res.x[]) > maxres)
                maxres = fabs (res.x[]);
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

            if (ibm[] <= 0.5)
                res.x[] = 0;

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

    restriction ({mu, rho, cm});
    struct Viscosity p = { mu, rho, dt };
    return mg_solve ((scalar *){u}, (scalar *){r},
                     residual_diffusion, relax_diffusion, &p, 
                     nrelax, res, minlevel = 1, 
                     tolerance = TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE);
}

