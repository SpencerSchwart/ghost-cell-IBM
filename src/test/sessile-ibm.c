#undef dv()
#define dv() (sq(Delta)*ibm[])

#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm.h"

#define LEVEL 6
#define MIN_LEVEL 3

double theta0;

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (0)

int main()
{
  size (2.);

  /**
  We shift the bottom boundary. */

  origin (0, - 0.26);
  init_grid (1 << LEVEL);
  
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */

  for (theta0 = 135; theta0 <= 165; theta0 += 15) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
}


event init (t = 0)
{

  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y;
  boundary ({phi});
  fractions (phi, ibm, ibmf);
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

event logfile (i++; t <= 20)
{

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;
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
    if (ibm[] < 1.)
      kappa[] = nodata;
  
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.4g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi);
}

#if 0
event adapt (i++)
{
  adapt_wavelet ({ibm,f,u}, (double[]){1.e-15, 1.e-8, 1.e-4, 1.e-4},
                 maxlevel = LEVEL, minlevel = MIN_LEVEL);
}
#endif
