#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../contact-ibm.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "view.h"

#define MAXLEVEL 7
#define MINLEVEL 3

#define R0 0.50
#define xc 0.
#define yc 0.575 

double t_end = 10;
double theta0;


u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

pf[bottom] = dirichlet(0);
p[bottom] = dirichlet(0);

int main() {
  //size (2.5);
  size(4.);

  /**
  We set the origin */

  origin (-L0/2., 0);

  init_grid (1 << MAXLEVEL);
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;
  TOLERANCE = 1e-5;
  /**
  We vary the contact_angle. */
#if 1
  const double angles[9] = {5,15,30,60,90,120,135,150,165};
  for (int i = 0; i < 9; ++i) {
    theta0 = angles[i];
  	const scalar c[] = angles[i]*pi/180.;
  	contact_angle = c;
    t_end = theta0 == 5? 60: theta0 == 15? 20: 15;
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
  We define the cylinder and the initial (half)-circular
  interface. */
  solid (ibm, ibmf, (sq(x - xc) + sq(y - yc) - sq(R0)));
  fractions_cleanup(ibm, ibmf, 1e-6);
  fraction (f, - (sq(x - xc) + sq(y - (yc + sqrt(2)/2)) - sq(R0)));
  fraction (ch, - (sq(x - xc) + sq(y - (yc + sqrt(2)/2)) - sq(R0)));
  foreach(reduction(+:v0))
    v0 += f[]*dv3();
}

event logfile (i++; t <= t_end)
{
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
 
 #if 0
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 5e-6)
    return true;
#endif

  double vreal = i == 0? v0: 0;
  foreach(reduction(+:vreal))
    vreal += f[]*sq(Delta);

  double verr = i == 0? 0: (vreal - v0)/v0 * 100;

  fprintf(stderr, "%d %g %g %g %g %g\n", i, t, theta0, v0, vreal, verr);
  fprintf(stdout, "%d %g %g %g %g %g\n", i, t, theta0, v0, vreal, verr);
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
  double R = s.volume/s.sum, V = statsf(f).sum;

  static FILE * fp2 = fopen("results", "w");

  fprintf (fp2, "%d %g %.5g %.5g %.3g %g %g\n", N, theta0, R, R/sqrt(V/pi), s.stddev, v0, V);

  fflush(fp2);
  if (theta0 == 165)
      fclose(fp2);
}

#if 0
event movie(i+=10,last){
    char name[80];
    sprintf(name, "movie-%g.mp4", theta0);
    view(fov=20, tx = 0, ty = -0.5);
    draw_vof ("f", lw=2);
    draw_vof("ibm", "ibmf",filled=-1);
    squares("f", linear = true, min = 0, max = 1);
    save(name);
}
#endif

#if 1
event adapt (i++) {
  scalar f1[], ibm1[];
  foreach() {
    f1[] = f[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1}, (double[]){1e-3, 1e-3}, 
    maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}
#endif

