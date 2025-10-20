#include "embed.h"
//#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../contact-embed.h"
#include "view.h"
#include "navier-stokes/perfs.h"

#define MAXLEVEL 8
#define MINLEVEL 5

double xc = 0;
double yc = 0;

double xd0 = 1.2;       // initial drop position
double ud0 = 1;         // initial drop velocity
double df = 1;          // drop diameter
double ds = 1.35;       // sphere diameter
double gravity = 0.039770;

double tmovie = 0.01;
double t_end = 1.5;
double theta0 = 95.1;

u.t[embed] = dirichlet(0);
u.n[embed] = dirichlet(0);

//p[top] = dirichlet(0);
//pf[top] = dirichlet(0);

u.n[left] = dirichlet(0);
u.t[left] = dirichlet(0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

f[bottom] = neumann(0);

int main() {
  size(4.);

  /**
  We set the origin */

  origin (-L0/2., 0);

  init_grid (1 << MAXLEVEL);
  
  rho1 = 1;
  rho2 = 0.001208;
  mu1 = 6.11e-4;
  mu2 = 1.09e-5;
  f.sigma = 0.059405;

  TOLERANCE = 1e-6;
  TOLERANCE_MU = 1e-3;
  
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
}

double v0 = 0;
event init (t = 0)
{
  /**
  We define the cylinder and the initial (half)-circular
  interface. */
  solid (cs, fs, (sq(x - xc) + sq(y - yc) - sq(ds*0.5)));
  fractions_cleanup(cs, fs, 1e-6);

#if 0
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);

  cm.refine = cm.prolongation = refine_cm_axi;
  cs.refine = cs.prolongation = fraction_refine;
  fm.x.refine = refine_face_x_axi;
  fm.y.refine = refine_face_y_axi;
  metric_embed_factor = axi_factor;

  restriction ({cs, fs, cm, fm});
#endif
  fraction (f, - (sq(x - xd0) + sq(y - yc) - sq(df*0.5)));

  foreach()
    u.x[] = -ud0*f[];
}

#if 0
event acceleration (i++) 
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= gravity;
}
#endif

event logfile (i++; t <= t_end)
{
  foreach_face()
    if (!fm.x[])
      uf.x[] = 0;

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  #if 0
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 5e-6)
    return true;
  #endif

  double vreal = i == 0? v0: 0;
  foreach(reduction(+:vreal))
    vreal += f[]*dv();

  if (i == 0) v0 = vreal;

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
  output_facets (f, fp);

 /**
  We compute the curvature only in full cells. */
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
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
event movie (t += tmovie) 
{
    char name[80];
    sprintf(name, "movie-%g.mp4", theta0);
    view(fov = 10, tx = -0.125, width = 2000, height = 2000);
    draw_vof ("f", lw = 2);
    draw_vof("cs", "fs",filled = -1);
    squares("f", min = 0, max = 1);
    save(name);
}
#endif

event interface (t = {0.2, 0.45, 0.9, 1.35})
{
    char name[80];
    sprintf(name, "%d-interface-%g.png", N, t);
    view(fov = 10, tx = -0.125, width = 2000, height = 2000);
    draw_vof ("f", lw = 2);
    draw_vof("cs", "fs",filled = -1);
    squares("f", min = 0, max = 1);
    save(name);

    sprintf (name, "%d-shape-%g", N, t);
    FILE * fp = fopen (name, "w");
    output_facets (f, fp);
    fclose (fp);
}

event adapt (i++)
{
  scalar f1[], cs1[];
  foreach() {
    f1[] = f[];
    cs1[] = cs[];
  }
  adapt_wavelet ({cs, f1}, (double[]){1e-3, 1e-3}, 
    maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}
