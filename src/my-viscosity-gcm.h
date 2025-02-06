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
            if (ibm[] > 0.5)
                u.x[] = (dt * a + r.x[] * sq(Delta)) /
                        (sq(Delta) * (rho[] + lambda.x) + avgmu);
            else if (ibm[] <= 0.5)
                u.x[] = 0;  //du = 0, since we dont want to change the imposed ghost cell velocity
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
            res.x[] = r.x[] - (rho[] + lambda.x) * u.x[] - dt * a /Delta;

            if (ibm[] <= 0.5)
                res.x[] = 0;

            if (fabs (res.x[]) > maxres)
                maxres = fabs (res.x[]);
        }
    }
    return maxres;
}

#undef lambda

double TOLERANCE_MU = 0.;

trace
mgstats viscosity (vector u, face vector mu, scalar rho, double dt,
                   int nrelax = 4, scalar * res = NULL)
{
    vector r[];

    foreach()
        foreach_dimension()
            r.x[] = rho[] * u.x[];

    restriction ({mu, rho});
    struct Viscosity p = { mu, rho, dt };
    return mg_solve ((scalar *){u}, (scalar *){r},
                     residual_diffusion, relax_diffusion, &p, 
                     nrelax, res, minlevel = 1, 
                     tolerance = TOLERANCE_MU ? TOLERANCE_MU : TOLERANCE);
}

