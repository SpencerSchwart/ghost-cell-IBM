#define CA 1
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm.h"

#define D0 1.

const int maxlevel = 8;
const int minlevel = 3;
const double L = 4.;

double t_end = 15.;
double theta0;

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

u.n[top] = neumann(0);
p[top]  = dirichlet(0);
pf[top] = dirichlet(0);

u.n[right] = neumann(0);
p[right]  = dirichlet(0);
pf[right] = dirichlet(0);

u.n[left] = neumann(0);
p[left]  = dirichlet(0);
pf[left] = dirichlet(0);

u.n[bottom] = neumann(0);
p[bottom]  = dirichlet(0);
pf[bottom] = dirichlet(0);

int main()
{
  size (L);
  origin (-3*L/4.,-L/4., 0);
  init_grid (1 << maxlevel);
  
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */
  TOLERANCE = 1e-5;

#if 0
  double angles[11] = {15,30,45.05,60,75,90,105,120,135.05,150,165};
  for (int i = 0; i < 11; i++) {
    theta0 = angles[i];
    t_end = theta0 == 15? 60: 15;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
#else
  t_end = 15;
  theta0 = 90;
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
#endif
}

double v0 = 0;

event init (t = 0)
{
  /**
  We define the inclined wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - x + 0.03);
  boundary ({phi});
  fractions (phi, ibm, ibmf);
  fraction (f, - (sq(x - 0) + sq(y - 0.) - sq(D0/2.)));
  fraction (ch, - (sq(x - 0) + sq(y - 0.) - sq(D0/2.)));

  v0 = real_volume(f);
}


event acceleration (i++)
{
    foreach_face(y)
      a.y[] += -0.1;

    foreach_face(x)
      a.x[] +=  0.1;
}

event logfile (i++; t <= t_end)
{

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;

  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += cr[]*sq(Delta);

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
  output_facets (ch, fp);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(cr).sum;

  static FILE * fp2 = fopen ("results", "w");

  fprintf (fp2, "%d %g %.5g %.3g %.4g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi, v0, V);

  fflush(fp2);
  if (theta0 == 165) // last case
      fclose (fp2);
}

event adapt(i++)
{
  scalar f1[], ibm1[];
  foreach() {
    f1[] = ch[];
    ibm1[] = ibm[];
  }
  adapt_wavelet({ibm1,f1, u}, (double[]){1e-3, 1e-3, 1e-4, 1e-4}, maxlevel = maxlevel, minlevel = minlevel);
}

