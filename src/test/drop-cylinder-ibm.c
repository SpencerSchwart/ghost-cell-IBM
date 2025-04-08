#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm.h"
#include "view.h"

#define R0 0.5
#define xc 0.
#define yc 0.575 

const double t_end = 15.;
double theta0;

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (0)


int main() {
  size (2.);

  /**
  We set the origin */

  origin (-1, 0);

  init_grid (64);
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */
#if 0
  for (theta0 = 30; theta0 <= 150; theta0 += 30) {
  	const scalar c[] = theta0*pi/180.;
  	contact_angle = c;
  	run();
  }
#else
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

  fraction (f, - (sq(x - xc) + sq(y - (yc+sqrt(2)/2)) - sq(R0)));
  v0 = real_volume(f);
}

event logfile (i++; t <= t_end)
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

  double vreal = real_volume(f), vreal2 = 0;
  foreach()
    vreal2 += cr[]*sq(Delta);

  fprintf(stderr, "%d %g %g %g %g %g\n", i, t, theta0, v0, vreal, vreal2);
}

#if 0
event volume (i++, t<=T){
  if (t==0) volume_vof_init = statsf (f).sum;

  char name[80];
  sprintf (name, "volume-mesh%d-angle%g.dat", N, theta0);
  static FILE * fp = fopen (name,"w");
  stats s = statsf (f);
  double erreur = ((volume_vof_init - s.sum)/volume_vof_init)*100;
  fprintf (fp, "%g %.5g %.5g\n", t, erreur, dt); 
}
#endif


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
  double R = s.volume/s.sum, V = real_volume(f);
  fprintf (stderr, "%d %g %.5g %.5g %.3g %g %g\n", N, theta0, R, R/sqrt(V/pi), s.stddev, v0, V);
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

