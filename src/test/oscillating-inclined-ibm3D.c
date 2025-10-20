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

#define L0 1.
#define D0 0.5

double Oh = 0.025;
double theta0=150;
double tend = 15;

#define MAXLEVEL 6
#define MINLEVEL 3

u.r[immersed] = neumann(0);
u.t[immersed] = neumann(0);
u.n[immersed] = dirichlet(0);

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

  /**
  We shift the bottom boundary. */

  origin (-L0/2.,-L0/2., 0);

  init_grid(1 << MAXLEVEL);
  
  TOLERANCE = 1e-5;
  TOLERANCE_MU = 1e-3;

  rho1 = rho2 = 1;
  f.sigma = 0.1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * D0);

#if 1
  double angles[2] = {50, 130};
  for (int i = 0; i < 2; i++) {
    theta0 = angles[i];
    tend = 12.001;
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
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  double minsize = L0/(1 << MAXLEVEL);

  vertex scalar phi[];
  foreach_vertex()
    phi[] = y - x > 0? 1: -1;
  boundary ({phi});
  fractions (phi, ibm, ibmf);
  boundary({ibm, ibmf});
  fraction (f, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
  fraction (ch, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));

  boundary({f});
  event("update_metric");

  foreach()
    foreach_dimension()
        u.x[] = 0;
  boundary({u});

  v0 = real_volume (f);
}

event logfile (i++; t <= tend)
{
  double v1 = i == 0? v0: 0;
  foreach(reduction(+:v1))
    v1 += cr[]*dv()/(ibm[] + SEPS);

  coord ns = {-0.5,0.5,0}; // normal of solid surface
  normalize(&ns);
  coord ts = {-ns.y, ns.x, 0};
  double hmax = -HUGE, rmax=-HUGE;
  foreach(reduction(max:hmax) reduction(max:rmax)) {
    if (ch[] > 1e-3 && ch[] < 1 && ibm[] && cr[] > 1e-3) {
        coord n = interface_normal (point, f), p;
        double alpha = plane_alpha (f[], n);
        plane_area_center (n, alpha, &p);
        coord cc = {x,y,z}, pp;
        foreach_dimension()
            pp.x = cc.x + p.x*Delta;

        double h = pp.x*ns.x + pp.y*ns.y + pp.z*ns.z;
        double r = pp.x*ts.x + pp.y*ts.y + pp.z*ts.z;
        if (h > hmax)
            hmax = h;
        if (r > rmax)
            rmax = r;
    }
  }

  double verr = i == 0? 0: (v1 - v0)/v0 * 100;

  fprintf(stdout, "%d %g %g %g %g %g %g %g\n", i, t, theta0, v0, v1, verr, hmax, rmax);
  fprintf(stderr, "%d %g %g %g %g %g %g %g\n", i, t, theta0, v0, v1, verr, hmax, rmax);

  if (i > 0) {
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-4)
    return true;
  }
}

#if 1
event movie (t += 0.1; t <= tend)
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

  view (camera="left", fov = 20, bg = {1,1,1}, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof ("ibm","ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof ("ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof ("ch", edges = true);
  squares ("ch");
  sprintf(name, "%g-chleft-%g.png", theta0, t);
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

event dump (t += 0.25)
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


#if 1
#if TREE
event adapt (i++) {
#if 1
  scalar f1[], ibm1[];
  foreach() {
    f1[] = f[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1}, (double[]){1e-3, 1e-3}, minlevel = MINLEVEL, maxlevel = MAXLEVEL);
#else
  adapt_wavelet ({f}, (double[]){1e-4}, minlevel = 3, maxlevel = MAXLEVEL);
#endif
}
#endif
#endif

/**
We compare $R/R_0$ to the analytical expression.

The accuracy is not as good (yet) as that of the [sessile
test](/src/test/sessile3D.c).

~~~gnuplot


reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 30,1 to 150,1 nohead dt 2
kappa(theta) = 2.*((2. + cos(theta))*(1. - cos(theta))**2/4.)**(1./3.)
R0(V) = (3.*V/(4.*pi))**(1./3.)
set xtics 30,30,150
plot 2./(kappa(x*pi/180.)) lw 3 lt -1 dt 2 lc rgb "black" t 'analytical', \
     'log' u 2:(2.*$3/R0($5)) w p pt 4 ps 3 lt -1 lw 3 t 'numerical'

~~~

## See also

* [3D Sessile drop on a domain boundary](/src/test/sessile3D.c)

*/



