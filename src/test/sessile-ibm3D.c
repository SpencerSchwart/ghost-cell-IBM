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
#include "profiling.h"
#include "view.h"

//#define T 1.5
#define L0 1.
#define D0 0.5

double t_end = 2;
double Oh = 0.05;
double theta0=150, volume_vof_init;

#define MAXLEVEL 6
#define MINLEVEL 3

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

u.n[left] = dirichlet(0);
u.n[back] = dirichlet(0);

f[left] = neumann(0);
f[back] = neumann(0);

p[left] = neumann(0);
p[back] = neumann(0);

pf[left] = neumann(0);
pf[back] = neumann(0);

int main()
{
  size (L0);
  init_grid(1 << MAXLEVEL);
  origin (0, -0.26, 0);

  TOLERANCE = 1e-5;
  TOLERANCE_MU = 1e-3;

#if 1
  mu1 = mu2 = 0.1;
  rho1 = rho2 = 1;
  f.sigma = 1.;
#else
  rho1 = rho2 = 1;
  f.sigma = 0.1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * D0);
#endif
#if 1
  double angles[7] = {15,30,60,90.05,120,150,165};
  for (int i = 0; i < 7; i++) {
    theta0 = angles[i];
    t_end = theta0 == 15? 5.001: 2.001;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
#else
    theta0 = 70;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
#endif
}

double v0;
event init (t = 0)
{
  if(!restore("dump-1.55")) {  
    /**
    We define the horizontal bottom wall and the initial (half)-circular
    interface. */
    
    vertex scalar phi[];
    foreach_vertex()
      phi[] = y + 0.00521;
    boundary ({phi});
    fractions (phi, ibm, ibmf);
    boundary({ibm, ibmf});
    restriction ({ibm, ibmf});
    fraction (f, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
    fraction (ch, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
    boundary({f,ch});
  
    foreach()
      foreach_dimension()
          u.x[] = 0;
    boundary({u});
  
    v0 = real_volume (f);
  }
  else {
    vertex scalar phi[];
    foreach_vertex()
      phi[] = y + 0.00521;
    boundary ({phi});
    fractions (phi, ibm, ibmf);
    boundary (all);
    v0 = 0.00815726; // for level 6
  }
}

event logfile (i++; t <= t_end)
{
  double v1 = 0;
  foreach(reduction(+:v1))
    v1 += cr[] * pow(Delta, dimension);

  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.) kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(cr).sum;

  if (statsf (kappa).stddev < 1e-2)
    return true;

  scalar pos[];
  position (ch, pos, {0,1,0});
  double hmax = statsf(pos).max;
  position (ch, pos, {1,0,0});
  foreach()
    if (ibm[] < 1.) pos[] = nodata;
  double rmax = statsf(pos).max;

  double verr = i == 0? 0: (v1 - v0)/v0 * 100;

  fprintf(stdout, "%d %g %g %g %g %g %g %d %g %g %g %g\n", 
      i, t, theta0, v0, v1, verr, avgitr, g_count, hmax, rmax, R, V);
  fprintf(stderr, "%d %g %g %g %g %g %g %d %g %g %g %g\n", 
      i, t, theta0, v0, v1, verr, avgitr, g_count, hmax, rmax, R, V);
  fflush(stdout);
}

#if DUMP
event dump (t += 0.5)
{
  scalar * dumplist = {u,p,f,ch,cr,ibm};
  char name[80];
  sprintf(name, "dump-%g", t);
  dump (file = name, list = dumplist);
}
#endif

#if 1
event movie (t += 0.10; t <= t_end)
{
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 30, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.299, tz = -2.748, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "ch", edges = true);
  mirror (n = {1,0,0}) {
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "ch", edges = true);
    }
  mirror (n = {0,0,1}) {
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "ch", edges = true);
    mirror (n = {1,0,0}) {
    cells (n = {0,1,0});
    draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
    draw_vof (c = "ch", fc = {0.647,0.114,0.176}, lw = 2);
    draw_vof (c= "ch", edges = true);
    }
  }
  char name[80];
  sprintf(name, "%g-movie-%g.png", theta0, t);
  save (name);
}
#endif


event end (t = end)
{
  /** At equilibrium we output the (almost constant) radius, volume, maximum velocity and time.*/
  FILE * fp = fopen ("results", "a");
  scalar kappa[];
  curvature (ch, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(cr).sum;
  fprintf (fp, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
     N, theta0, R, s.stddev, V, normf(u.x).max, t);
  fflush (fp);
  fclose (fp);
}

event adapt (i++) {
  scalar f1[], ibm1[];
  foreach() {
    f1[] = f[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1}, (double[]){1e-3, 1e-3}, 
    maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

