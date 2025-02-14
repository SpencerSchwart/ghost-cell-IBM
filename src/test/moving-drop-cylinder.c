/**
# droplet on an embedded cylinder


~~~gnuplot Equilibrium shapes for $30^\circ \leq \theta \leq 150^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1:1]
set yrange [0:]
f0(x) = sqrt(0.5**2 - x**2) + 0.575
f(x)  = sqrt(0.5**2 - x**2) - 0.575

plot 'out' w l lt -1 lw 3 lc rgb "blue" t 'Numerical solution',\
  f0(x) with filledcurves above y1 = 0.575 fc "black" t 'cylinder',\
 -1*f(x) with filledcurves below y1 = 0.575 fc "black" t ''
 
set term pop
~~~
*/

#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../my-two-phase.h"
//#include "navier-stokes/conserving.h"
#include "../my-tension.h"
#include "../contact-ibm.h"
#include "view.h"

#define R0 0.5
#define L0 4.
#define xc 2.
#define yc 2.
#define T 5

double theta0, volume_vof_init;
int LEVEL = 8;
int MIN_LEVEL = 4;
const double A = 0.1 * R0;
const double freq = 0.1;
const double t_start = T;
const double t_end = T + 15;

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (t < t_start? 0: A*2*M_PI*freq*sin(2*M_PI*freq*(t - t_start)))

int main() {
  size (L0);

  /**
  We set the origin */

  //origin (-L0/2, -L0/2);

  init_grid (1 << LEVEL);
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */

  theta0 = 60;
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
}


event init (t = 0)
{
  /**
  We define the cylinder and the initial (half)-circular
  interface. */
  double ypos = t < t_start? -A :-A*cos(2*M_PI*freq*(t - t_start));
  solid (ibm, ibmf, sq(x - xc) + sq(y - yc - ypos) - sq(R0));

  fraction (f, - (sq(x - xc) + sq(y - (yc+sqrt(2)/2)) - sq(R0)));
}


event moving_cylinder (i++)
{
  double ypos = t < t_start? -A :-A*cos(2*M_PI*freq*(t - t_start));
  solid (ibm, ibmf, sq(x - xc) + sq(y - yc - ypos) - sq(R0));
}

event logfile (i++; t <= t_end)
{
    double lvolume = 0;
    foreach()
        lvolume += f[] * ibm[] * Delta;
    fprintf (stderr, "%d %g %d %d %d %d %d %d %g\n",
             i, t, mgpf.i, mgpf.nrelax, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, lvolume);
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

  output_facets (f, stdout);

 /**
  We compute the curvature only in full cells. */
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.5g %.3g\n", N, theta0, R, R/sqrt(V/pi), s.stddev);
}

event movie(i+=10,last){
  view(fov=20, tx = -0.5, ty = -0.5);
  draw_vof ("f", lw=2);
  draw_vof("ibm", "ibmf",filled=-1);
  squares("p", linear = true, min = 0, max = 1);
  save("movie.mp4");
}

event adapt (i++) {
  adapt_wavelet ({ibm,f,u}, (double[]){1.e-5,1e-2,3e-4,3e-4},
		 maxlevel = LEVEL, minlevel = MIN_LEVEL);
}

/**
![Relaxation toward a $120^\circ$ contact angle.](droplet-cylinder-embed/movie.mp4)
*/

// fixme: Comparison to theory is missing (add soon) 
