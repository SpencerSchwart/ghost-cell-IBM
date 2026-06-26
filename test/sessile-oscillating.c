/**
A 2D sessile droplet oscillating on a flat, horizontal plate until it reaches an
equilibirum position, which is based on the contact angle. The oscillation occurs
from low Ohnesorge number and initalizing the droplet at 90 degrees instead of
its equilibirum contact angle (60 or 120).

This file uses Basilisk's built in contact.h to impose the boundary condition
and is used as a reference. */

#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "vof.h"

/** 
Simulation variables */

int maxlevel = 6;
int minlevel = 4;

double df = 1;         // droplet diameter
double theta0 = 0;     // contact angle
double Oh = 0.005;     // Ohnesorge number
double v0 = 0;         // initial volume

double t_end = 40;

FILE * outf; // out file

/**
Boundary conditions. */

vector h[];
h.t[bottom] = contact_angle (theta0*pi/180.);

u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);

u.n[top] = dirichlet(0);
p[top]   = dirichlet(0);
pf[top]  = dirichlet(0);

int main(int argc, char* argv[])
{
  size (2.);

  if (argc > 1)
    theta0 = atof(argv[1]);
  if (argc > 2)
    t_end = atof(argv[2]);
  if (argc > 3)
    maxlevel = atoi(argv[3]);

  origin (0, 0);
  init_grid (1 << maxlevel);
  
  TOLERANCE = 1e-3;

  /**
  We use a constant fluid properties. */
  rho1 = rho2 = 1;
  f.sigma = 1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

  f.height = h;

  /**
  We vary the contact_angle if theta0 != 0. */
  if (theta0) 
    run();
  else {
    double angles[2] = {110,70};
    for (int i = 0; i < 2; i++) {
      theta0 = angles[i];
      run();
    }
  }
}


event init (t = 0)
{
  /**
  Initalize the droplet and calculate the initial volume. */
  fraction (f, - (sq(x) + sq(y) - sq(df/2.)));

  foreach()
    v0 += f[]*dv();
}


event logfile (i++; t <= t_end)
{
  /**
  Calculate volume, maximum capillary number, and maximum droplet height. */

  double vreal = 0, ca = 0, hf = 0;
  foreach(reduction(+:vreal) reduction(max:ca) reduction(max:hf)) {

    vreal += f[]*sq(Delta);

    double umag = sqrt(sq(u.x[]) + sq(u.y[]));
    if (umag > ca)
      ca = umag;

    if (f[] > 0 && f[] < 1) {
      coord n = interface_normal(point, f), mp;
      double alpha = plane_alpha(f[], n);
      plane_area_center (n, alpha, &mp);

      mp.y = y + Delta*mp.y;

      if (mp.y > hf)
        hf = mp.y;
    }
  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  /**
  Get the contact line position. */

  double xcl = 0;
  foreach_boundary(bottom) {

    if (f[] > 1e-6 && f[] < 1-1e-6) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);

      coord pint[2];
      int count = facets (n, alpha, pint);

      /**
      Find the point which intersects the bottom cells face/boundary (y = -0.5). */

      double xcl0 = 0;
      for (int i = 0; i < count; i++) 
        if (pint[i].y <= -0.5 + 1e-6)
            xcl0 = x + pint[i].x*Delta;

      if (xcl0 > xcl)
        xcl = xcl0;
    }
  }

  if (i == 0) {
    char name[80];
    sprintf(name, "outs/out-%g-%d", theta0, maxlevel);
    outf = fopen(name, "w");
  }

  fprintf(outf, "%d %g %g %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verror, xcl, ca, hf);
  fflush(outf);
}


/**
Given the radius of curvature $R$ and the volume of the droplet $V$,
this function returns the equivalent contact angle $\theta$ verifying
the equilibrium solution
$$
V = R^2(\theta - \sin\theta\cos\theta)
$$
*/

double equivalent_contact_angle (double R, double V)
{
  double x0 = 0., x1 = pi;
  while (x1 - x0 > 1e-4) {
    double x = (x1 + x0)/2.;
    double f = V - sq(R)*(x - sin(x)*cos(x));
    if (f > 0.)
      x0 = x;
    else
      x1 = x;
  }
  return (x0 + x1)/2.;
}


event end (t = end)
{
  fclose(outf);

  scalar kappa[];
  curvature (f, kappa);

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;

  char name[80];
  sprintf(name, "results/results-%g-%d", theta0, maxlevel);
  FILE * fp2 = fopen (name, "w");

  fprintf (fp2, "%d %g %.5g %.3g %.4g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi, v0, V);
  fclose(fp2);
}

/**
We trick the adapt function to keep the liquid interface at the maximum level
of refinement. */
event adapt (i++)
{
  scalar f1[];
  foreach()
    f1[] = f[];

  adapt_wavelet ({f1, u}, (double[]){1e-3, 1e-2, 1e-2}, 
    minlevel = minlevel, maxlevel = maxlevel);
}

