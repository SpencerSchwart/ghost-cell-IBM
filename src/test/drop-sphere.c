#define CA 1
#define DUMP 1

#include "grid/octree.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm.h"
#include "navier-stokes/perfs.h"
#include "view.h"

#if TRACE
#include "profiling.h"
#endif

#define MAXLEVEL 7
#define MINLEVEL 3

double tend = 15;
#define L0 2.1
#define D 1.

coord ci = {0,0,0};
coord di = {0,sqrt(2.)/2.,0}; // sqrt(2)/2

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

double theta0;

int main() {
  size(L0);
  init_grid (1 << MINLEVEL);
  origin (0, -0.575, 0);
  mu1 = mu2 = 0.1;
  f.sigma = 1.;

  TOLERANCE = 1.e-5;

#if 1
  //const double theta0_array[9] = {5,15,30,60,90,120,135,150,165};
  const double theta0_array[8] = {15,30,60,90,120,135,150,165};
  for (int i = 0; i < 8; ++i) { // CHANGE to 8 or 9
    theta0 = theta0_array[i];
  	const scalar c[] = theta0_array[i]*pi/180.;
  	contact_angle = c;
    if (theta0 == 15)
        tend = 25.001;
    else
        tend = 15.001;
  	run();
  }
#else
    theta0 = 15;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
#endif
}

double v0 = 0;
event init (t = 0) 
{
  scalar ibmFake[], fFake[];
  face vector ibmfFake[];

  int count = 0;
  astats st;
  do {
    solid (ibmFake, ibmfFake, sq(x - ci.x) + sq(y - ci.y) + sq(z - ci.z) - sq(D/2));
    fraction (fFake, - (sq(x - di.x) + sq(y - di.y) + sq(z - di.z) - sq(D/2)));
    st = adapt_wavelet ({ibmFake, fFake}, (double[]){1e-5, 1e-5}, MAXLEVEL, MINLEVEL);
    count++;
  } while ((st.nc || st.nf) && count < 30);

  solid (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) + sq(z - ci.z) - sq(D/2));
  fraction (f, - (sq(x - di.x) + sq(y - di.y) + sq(z - di.z) - sq(D/2)));
  fraction (ch, - (sq(x - di.x) + sq(y - di.y) + sq(z - di.z) - sq(D/2)));

  v0 = real_volume(f);
}

event logfile (i++; t <= tend)
{
  double v1 = i == 0? v0: 0;
  foreach(reduction(+:v1))
    v1 += cr[]*dv()/(ibm[] + SEPS);

  double verr = i == 0? 0: (v1 - v0)/v0 * 100;
  fprintf(stdout, "%d %g %g %g %g %g\n", i, t, theta0, v0, v1, verr);
  fprintf(stderr, "%d %g %g %g %g %g\n", i, t, theta0, v0, v1, verr);

  /**
  If the standard deviation of curvature falls below $10^{-2}$, we stop the computation
  (convergence has been reached). */  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-2)
    return true;
}

#if DUMP
event dumps (t += 2.5; t <= tend)
{
  scalar * dumplist = {u,p,f,ibm};
  char name[80];
  sprintf(name, "%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}
#endif

#if 1
event movie (t += 0.1; t <= tend)
{
  char name[80];
  double fov = 18, ty = -0.2;

  view (camera="left", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("f", edges = true);
  sprintf(name, "%g-left-%g.png", theta0, t);
  save (name);

  view (camera="left", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("cr", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("cr", edges = true);
  sprintf(name, "%g-crleft-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells ();
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("cr", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("cr", edges = true);
  squares ("cr");
  sprintf(name, "%g-crright-%g.png", theta0, t);
  save (name);

  view (camera="left", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  sprintf(name, "%g-chleft-%g.png", theta0, t);
  save (name);

  view (camera="back", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("f", edges = true);
  squares ("f");
  sprintf(name, "%g-back-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells ();
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("f", edges = true);
  squares ("f");
  sprintf(name, "%g-right-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells ();
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  squares ("ch");
  sprintf(name, "%g-chright-%g.png", theta0, t);
  save (name);

  view (camera="back", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  squares ("ch");
  sprintf(name, "%g-chback-%g.png", theta0, t);
  save (name);
}
#endif

event end (t = end)
{
  /** At equilibrium we output the (almost constant) radius, volume, maximum velocity and time.*/
  FILE * fp = fopen ("results", "a");
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
        kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(cr).sum;

  fprintf (fp, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
     N, theta0, R, s.stddev, V, normf(u.x).max, t);
  fflush(fp);
  fclose(fp);

  scalar * dumplist = {u,p,f,ibm};
  char name[80];
  sprintf(name, "%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}

event adapt (i++) 
{
  scalar f1[], ibm1[];
  foreach() {
    f1[] = f[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1}, (double[]){1e-3, 1e-3}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
}
