/**
A 2D sessile droplet spreading on a flat, horizontal plate until it reaches an
equilibirum position, which is based on the contact angle. */

#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

/** 
Simulation variables */

const int maxlevel = 6;
const int minlevel = 4;

double df = 1;         // droplet diameter
double theta0 = 0;     // contact angle
double Oh = 0.04;      // Ohnesorge number
double v0 = 0;         // initial volume

double t_end = 20;

FILE * outf; // out file

/** 
Boundary conditions */

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

u.n[top] = neumann(0);
p[top] = dirichlet(0);
pf[top] = dirichlet(0);

int main(int argc, char* argv[])
{
  size (2.);

  /**
  We shift the bottom boundary. */
  origin (0, -0.254);
  init_grid (1 << maxlevel);
  
  TOLERANCE = 1e-3;

  if (argc > 1)
    theta0 = atof(argv[1]);
  if (argc > 2)
    t_end = atof(argv[2]);
  if (argc > 3)
    Oh = atof(argv[3]);

  /**
  We use a constant fluid properties. */
  rho1 = rho2 = 1;
  f.sigma = 1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

  /**
  We vary the contact_angle if theta0 != 0. */
  if (theta0) {
    f.wetting.theta_s = theta0;
    run();
  } 
  else {
    double angles[11] = {15,30,45.05,60,75,90.05,105,120,135.05,150,165};
    double tends[11] =  {90,30,20,20,20,15,15,15,10,10,10};
    for (int i = 0; i < 11; i++) {
      theta0 = angles[i];
      t_end = tends[i];

      if (theta0 >= 150)
        Oh = 0.16;

      mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

      f.wetting.theta_s = theta0;
      run();
    }
  }
}

void solid_domain (scalar cs, face vector fs)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y;
  boundary ({phi});
  fractions (phi, cs, fs);
}

event init (t = 0)
{
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  solid_domain(cs, fs);
  fraction (f, - (sq(x) + sq(y) - sq(df/2.)));
  fraction (ch, - (sq(x) + sq(y) - sq(df/2.)));

  v0 = real_volume(f, cs, fs);
}

event logfile (i++; t <= t_end)
{

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1. || f[] > cs[] - 1e-6)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-5)
    return true;

  /**
  Calculate volume, contact line position, and maximum capillary number. */

  double vreal = 0, xcl = 0, ca = 0;
  foreach() {

    vreal += f[]*sq(Delta);

    if (on_interface(cs) && f[] > INT_TOL && f[] < cs[] - INT_TOL) {
      coord pint;
      if (is_contact_cell_bc(point, f, cs, fs, pint = &pint)) {
        pint = local_to_global_coord(point, pint);
        if (pint.x > xcl)
          xcl = pint.x;
      }
    }

    if (gc[] && cs[] == 1 && f[]) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]));
      if (umag > ca)
        ca = umag;
    }
  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  if (i == 0) {
    char name[80];
    sprintf(name, "outs/out-%g", theta0);
    outf = fopen(name, "w");
  }

  fprintf(outf, "%d %g %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verror, xcl, ca);
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

  /**
  At the end, we output the equilibrium shape. */
  
  char name[80];
  sprintf (name, "shapes/shape-%g", theta0);
  FILE * fp = fopen (name, "w");
  output_facets_contact (f, ch, cs, fp);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1. || f[] > cs[] - 1e-6)
      kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;

  sprintf(name, "results/results-%g", theta0);
  FILE * fp2 = fopen (name, "w");

  fprintf (fp2, "%d %g %.5g %.3g %.4g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi, v0, V);
  fflush(fp2);
  fclose(fp2);
}

/**
We use an intermediate field f1 to trick the adapt function into always
keeping the liquid interface at the maximum resolution. */

#if 1
event adapt (i++)
{
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1, u,}, (double[]){1e-3, 1e-2, 1e-2},
                           maxlevel = maxlevel, minlevel = minlevel);
}
#endif
