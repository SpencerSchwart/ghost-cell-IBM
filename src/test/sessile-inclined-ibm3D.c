#define CA 1
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

#define L0 1.5

#define R0 0.25

double t_end = 2.;
double theta0=150;

#define MAXLEVEL 7
#define MINLEVEL 3

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

u.n[left] = dirichlet(0);
u.n[back] = dirichlet(0);

p[back]  = neumann(0);
pf[back] = neumann(0);

p[front]  = dirichlet(0);
pf[front] = dirichlet(0);

u.n[front] = neumann(0);
u.t[front] = neumann(0);

p[top]  = dirichlet(0);
pf[top] = dirichlet(0);

u.n[top] = neumann(0);
u.t[top] = neumann(0);

int main()
{
  size (L0);

  origin (-L0/2.,-L0/2., 0);

  init_grid(1 << MAXLEVEL);
  
  TOLERANCE = 1e-5;
  TOLERANCE_MU = 1e-3;

  // We use a constant viscosity. 
  mu1 = mu2 = 0.1;
  rho1 = rho2 = 1.;
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */
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
    theta0 = 15;
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
#endif
}

double v0;
event init (t = 0)
{
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
  phi[] = y - x + 0.12;
  boundary ({phi});
  fractions (phi, ibm, ibmf);
  boundary({ibm, ibmf});
  fraction (f, - (sq(x - 0.06/sqrt(2)) + sq(y + 0.06/sqrt(2)) + sq(z) - sq(R0)));
  fraction (ch, - (sq(x - 0.06/sqrt(2)) + sq(y + 0.06/sqrt(2)) + sq(z) - sq(R0)));

  boundary({f, ch});

  foreach()
    foreach_dimension()
        u.x[] = 0;
  boundary({u});

  v0 = real_volume (f);
}

event logfile (i++; t <= t_end)
{
  double v1 = i == 0? v0: 0;
  foreach(reduction(+:v1))
    v1 += cr[]*dv();

  double verr = i == 0? 0: (v1 - v0)/v0 * 100;
  fprintf(stdout, "%d %g %g %g %g %g\n", i, t, theta0, v0, v1, verr);
  fprintf(stderr, "%d %g %g %g %g %g\n", i, t, theta0, v0, v1, verr);

  if (i > 0) {
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-2)
    return true;
  }
}

#if 1
event movie (t += 0.05; t <= t_end)
{
  char name[80];

  view (camera="front", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("cr", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("cr", edges = true);
  sprintf(name, "%g-crfront-%g.png", theta0, t);
  save (name);

  view (camera="back", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("cr", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("cr", edges = true);
  squares ("cr");
  sprintf(name, "%g-crback-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("cr", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("cr", edges = true);
  squares ("cr");
  sprintf(name, "%g-crright-%g.png", theta0, t);
  save (name);

  view (camera="front", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  sprintf(name, "%g-chfront-%g.png", theta0, t);
  save (name);

  view (camera="back", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  squares ("ch");
  sprintf(name, "%g-chback-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  squares ("ch");
  sprintf(name, "%g-chright-%g.png", theta0, t);
  save (name);

  view (camera="front", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("f", edges = true);
  sprintf(name, "%g-ffront-%g.png", theta0, t);
  save (name);

  view (camera="back", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("f", edges = true);
  squares ("f");
  sprintf(name, "%g-fback-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("f", edges = true);
  squares ("f");
  sprintf(name, "%g-fright-%g.png", theta0, t);
  save (name);

}
#endif

event dump (t += 0.5)
{
  scalar * dumplist = {u,p,f,ibm,cr,ch};
  char name[80];
  sprintf(name, "dumps/%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}

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
}

event adapt (i++) {
  scalar f1[], ibm1[];
  foreach() {
    f1[] = f[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1}, (double[]){1e-1, 1e-15}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
}

