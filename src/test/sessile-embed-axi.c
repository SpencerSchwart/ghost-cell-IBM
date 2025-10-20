#define CA 0

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "../my-two-phase.h"
#include "tension.h"
#include "../contact-embed.h"

#define D0 1.

const int maxlevel = 6;

double t_end = 15;
double theta0, Oh = 0.05;

u.t[embed] = dirichlet(0);
u.n[embed] = dirichlet(0);

u.n[bottom] = dirichlet(0);

int main()
{
  size (2.);

  /**
  We shift the bottom boundary. */

  origin (-0.26, 0);
  init_grid (1 << maxlevel);
  
  TOLERANCE = 1e-5;

  /**
  We use a constant viscosity. */
#if 1
  mu1 = mu2 = 0.1;
  f.sigma = 1.;
#else
  rho1 = rho2 = 1;
  f.sigma = 0.1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * D0);

#endif
  /**
  We vary the contact_angle. */

#if 1
  double angles[11] = {15,30,45,60,75,90.05,105,120,135,150,165};
  for (int i = 0; i < 11; i++) {
    theta0 = angles[i];
    t_end = theta0 == 15? 50: 15;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
#else
    t_end = 40;
    theta0 = 165;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
#endif
}

double v0 = 0;

event init (t = 0)
{

  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = x;
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - (sq(x) + sq(y) - sq(D0/2.)));

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);

  cm.refine = cm.prolongation = refine_cm_axi;
  cs.refine = cs.prolongation = fraction_refine;
  fm.x.refine = refine_face_x_axi;
  fm.y.refine = refine_face_y_axi;
  metric_embed_factor = axi_factor;

  restriction ({cs, fs, cm, fm});

  v0 = 0;
  foreach()
    v0 += cs[]? f[]*sq(Delta)*cm[]: 0;
}

face vector fm0[];
scalar cm0[];

event logfile (i++; t <= t_end)
{

  foreach_face()
    fm0.x[] = fm.x[];
  foreach()
    cm0[] = cm[];

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;

  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += cs[]? f[]*sq(Delta)*cm[]: 0;

#if 0
  scalar pos[];
  position (cr, pos, {0,1});
  double hmax = statsf(pos).max;
#else
  double hmax = 0;
#endif

  double perror = i == 0? 0: (vreal - v0)/v0 * 100;

  fprintf(stderr, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, perror, hmax);
  fprintf(stdout, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, perror, hmax);
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
  /**
  At the end, we output the equilibrium shape. */
  
  char name[80];
  sprintf (name, "shape-%g", theta0);
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = pi*sq(D0/.2)*statsf(f).sum;

  static FILE * fp2 = fopen ("results", "w");

  fprintf (fp2, "%d %g %.5g %.3g %.4g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi, v0, V);

  fflush(fp2);
  if (theta0 == 165) // last case
      fclose (fp2);
}

