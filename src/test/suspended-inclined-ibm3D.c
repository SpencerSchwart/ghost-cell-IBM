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

#define D0 1

u.t[immersed] = neumann(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = neumann(0);

#if 0
u.n[left] = dirichlet(0);
u.t[left] = neumann(0);
u.r[left] = neumann(0);

u.n[back] = dirichlet(0);
u.t[back] = neumann(0);
u.r[left] = neumann(0);

f[left] = neumann(0);
f[back] = neumann(0);

p[left] = neumann(0);
p[back] = neumann(0);

pf[left] = neumann(0);
pf[back] = neumann(0);

u.n[top] = neumann(0);
p[top]  = dirichlet(0);
pf[top] = dirichlet(0);
#else

u.n[top] = neumann(0);
p  [top] = dirichlet(0);
pf [top] = dirichlet(0);

u.n[bottom] = neumann(0);
p  [bottom] = dirichlet(0);
pf [bottom] = dirichlet(0);

u.n[right] = neumann(0);
p  [right] = dirichlet(0);
pf [right] = dirichlet(0);

u.n[left] = neumann(0);
p  [left] = dirichlet(0);
pf [left] = dirichlet(0);

u.n[back] = dirichlet(0);
u.t[back] = neumann(0);
f  [back] = neumann(0);

#endif


int max_level = 9;
int min_level = 3;
double L = 4.;
double t_end = 20;
double amp = 0.e-2;
double dia = 5.e-3;
double gpert = 4.9;

double femax = 0.01;
double uemax = 0.001;

double theta0 = 90.01;
double rhol   = 1000;
double rhog   = 1.2;
double mul    = 1.e-3;
double mug    = 1.8e-5;
double sigma  = 0.07;

double omegac = 0.5367;
double time_pert = 1.0434;


int main()
{
  size (L);

  //origin (-L/2.,-L/2., 0);
  origin (-(4.5)*L/(5.5),-L/(5.5), 0); 

  init_grid(1 << min_level);
  
  TOLERANCE = 1e-5;

  rho1 = rhol/rhol;
  rho2 = rhog/rhol;
  mu1 = mul/sqrt(rhol*dia*sigma);
  mu2 = mug/sqrt(rhol*dia*sigma);

  f.sigma = sigma/sigma;

  const scalar c[] = theta0*pi/180.;
  contact_angle = c;

  run();
}

void solid_plate (scalar s, face vector sf)
{
  vertex scalar phi[];
  foreach_vertex() 
    phi[] = max(y - x, z - 1);
  boundary ({phi});
  fractions (phi, s, sf);
  boundary({s, sf});
}

double v0;
event init (t = 0)
{
if (!restore("stampede/dump-0.4")) {
  scalar ibm1[], f1[];
  face vector ibmf1[];

  int count = 0;
  astats as;
  do {
    solid_plate (ibm1, ibmf1);
    fraction (f1, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
    as = adapt_wavelet ({ibm1, f1}, (double[]){1e-5,1e-5}, minlevel = min_level, maxlevel = max_level);
    count++;
  } while ((as.nc || as.nf) && count < 50);

  solid_plate (ibm, ibmf);
  fraction (f, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
  fraction (ch, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));

  boundary({f, ch});
  event("update_metric");

  foreach()
    foreach_dimension()
        u.x[] = 0;
  boundary({u});

  v0 = real_volume (f);
}
else {
   solid_plate (ibm, ibmf);   
   boundary(all);
   v0 = real_volume (f);
}
}

event acceleration (i++)
{
  if (t <= time_pert) {
    foreach_face(y)
      a.y[] += (t > time_pert? 0:
               -sqrt(2)/2.0*gpert*rhol*sq(dia)/sigma);

    foreach_face(x)
      a.x[] += (t > time_pert? 0:
               sqrt(2)/2.0*gpert*rhol*sq(dia)/sigma);
  }
}

event logfile (i++; t <= t_end)
{

  foreach_face()
    if (!fm.x[])
      uf.x[] = 0;

  double v1 = i == 0? v0: 0;
  foreach(reduction(+:v1))
    v1 += cr[]*dv();

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

  if (i == 0) { hmax = D0/2., rmax = D0/2.; }

  double verr = i == 0? 0: (v1 - v0)/v0 * 100;

  fprintf(stdout, "%d %g %g %g %g %g %g %g\n", i, t, theta0, v0, v1, verr, hmax, rmax);
  fprintf(stderr, "%d %g %g %g %g %g %g %g\n", i, t, theta0, v0, v1, verr, hmax, rmax);
}

#if 1
event movie (t += 0.1; t <= t_end)
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

event dump (t += 0.1)
{
  scalar * dumplist = {u,p,f,ibm,cr,ch};
  char name[80];
  sprintf(name, "dump-%g", t);
  dump (file = name, list = dumplist);
}

#if 1
event adapt (i++) {
  scalar ch1[], ibm1[];
  foreach() {
    ch1[] = f[]; 
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, ch1, u}, (double[]){1e-3, 1e-3, 3e-2, 3e-2, 3e-2}, 
    minlevel = min_level, maxlevel = max_level);
  foreach_face()
    if (!fm.x[])
      uf.x[] = 0;
}
#endif

