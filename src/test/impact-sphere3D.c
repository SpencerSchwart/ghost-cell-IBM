#define CA 1
#include "grid/octree.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../contact-ibm.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#if TRACE
#include "profiling.h"
#endif

#define MAXLEVEL 8
#define MINLEVEL 5

double xc = 0;
double yc = 0;
double zc = 0;

double xd0 = 1.9;       // initial drop position
double ud0 = 1;         // initial drop velocity
double df = 1;          // drop diameter
double ds = 2.7;       // sphere diameter
double gravity = 0.02308;

double tmovie = 0.0125;
double tdump = 0.1;
double t_end = 12;
double theta0 = 30.05;

u.n[immersed] = dirichlet(0);
u.t[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

p  [right] = dirichlet(0);
pf [right] = dirichlet(0);
u.n[right] = neumann(0);

p  [left] = dirichlet(0);
pf [left] = dirichlet(0);
u.n[left] = neumann(0);

int main() {
  size(10);

  /**
  We set the origin */

  origin (-ds/3, 0);

  init_grid (1 << MINLEVEL);
  
  rho1 = 1;
  rho2 = 0.001532;
  mu1 = 1.28E-3;
  mu2 = 1.07E-5;
  f.sigma = 0.012987;

  TOLERANCE = 1e-5;
  
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
}

event init (t = 0)
{
  /**
  We define the cylinder and the initial (half)-circular
  interface. */
  scalar ibm1[], f1[];
  face vector ibmf1[];

  int count = 0;
  astats st;
  do {
    solid (ibm1, ibmf1, (sq(x - xc) + sq(y - yc) + sq(z - zc) - sq(ds*0.5)));
    fraction (f1, - (sq(x - xd0) + sq(y - yc) + sq(z - zc) - sq(df*0.5)));
    st = adapt_wavelet ({ibm1, f1}, (double[]){1e-5, 1e-5}, MAXLEVEL, MINLEVEL);
    count++;
  } while ((st.nc || st.nf) && count < 40);

  solid (ibm, ibmf, (sq(x - xc) + sq(y - yc) + sq(z - zc) - sq(ds*0.5)));
  fractions_cleanup(ibm, ibmf, 1e-4);

  fraction(f,  - (sq(x - xd0) + sq(y - yc) + sq(z - zc) - sq(df*0.5)));
  fraction(ch, - (sq(x - xd0) + sq(y - yc) + sq(z - zc) - sq(df*0.5)));
  fraction(cr, - (sq(x - xd0) + sq(y - yc) + sq(z - zc) - sq(df*0.5)));

  foreach()
    u.x[] = -ud0*f[];
}

#if 1
event acceleration (i++) 
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= gravity;
}
#endif

double v0 = 0;
event logfile (i++; t <= t_end)
{
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += cr[]*sq(Delta)*cm[];

  if (i == 1) v0 = vreal;

  double hf = 0;
  foreach_boundary(bottom, reduction(max:hf)) {
    if (f[] > 1e-6 && f[] < 1-1e-6 && ibm[] == 1) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double htemp = x + p.x*Delta - 0.5*ds;
      if (htemp > hf)
        hf = htemp;
    }
  }

  double verr = i == 0? 0: (vreal - v0)/v0 * 100;

  fprintf(stderr, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verr, hf);
  fprintf(stdout, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verr, hf);
}

event dump (t += tdump)
{
  scalar * dumplist = {u,p,f,ibm,cr,ch};
  char name[80];
  sprintf(name, "dump-%g", t);
  dump (file = name, list = dumplist);
}

#if 1
event movie (t += tmovie) 
{
    char name[80];
    sprintf(name, "%d-movie-%g.png", N, t);
    view(camera="front", fov = 11, tx = -0.075, width = 1000, height = 1000);
    draw_vof ("ch", lw = 2, fc={1, 0.522, 0.439});
    draw_vof("ibm", "ibmf",filled = -1);
    squares("ch", min = 0, max = 1);
    save(name);
}
#endif

#if 1
event adapt (i++)
{
  scalar f1[], ibm1[];
  foreach() {
    f1[] = ch[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1, u}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-3}, 
    maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}
#endif
