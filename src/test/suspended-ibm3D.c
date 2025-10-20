#define CA 1
#include "grid/octree.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../my-two-phase.h"
#include "../my-conserving.h"
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

int max_level = 8;
int min_level = 5;
double L = 2.;
double t_end = 20;
double amp = 0.e-2;
double dia = 5.e-3;
double gpert = 4.9;

double femax = 0.01;
double uemax = 0.001;

double theta0 = 90;
double rhol   = 1000;
double rhog   = 1.2;
double mul    = 1.e-3;
double mug    = 1.8e-5;
double sigma  = 0.07;

double omegac    = 0.5367;
double time_pert = 1.0434;

int main()
{
  size (L);

  origin (0, -0.26, 0);
  
  init_grid(1 << min_level);
  
  TOLERANCE = 1e-5;
  TOLERANCE_MU = 1e-5;

  rho1 = rhol/rhol;
  rho2 = rhog/rhol;
  mu1 = mul/sqrt(rhol*dia*sigma);
  mu2 = mug/sqrt(rhol*dia*sigma);

  f.sigma = sigma/sigma;

  run();
}

void solid_domain (scalar s, face vector sf)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y;
  boundary ({phi});
  fractions (phi, s, sf);
  boundary({s, sf});
}

double v0;
event init (t = 0)
{
  if(!restore("70-dump-1.55")) {  
    
    scalar ibm1[], f1[];
    face vector ibmf1[];

    astats st; int count = 0;
    do {
      solid_domain(ibm1, ibmf1);
      fraction (f1, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
      st = adapt_wavelet ({ibm1, f1}, (double[]){1e-5, 1e-5}, max_level, min_level);
      count++;
    } while ((st.nc || st.nf) && count < 50);

    solid_domain (ibm, ibmf);
    fraction (f, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
    fraction (ch, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
    foreach()
      foreach_dimension()
          u.x[] = 0;
    boundary({u, ibm, f});
  
    v0 = real_volume (f);
  }
  else {
    solid_domain (ibm, ibmf);
    boundary (all);
    v0 = 0.00815726; // for level 6?
  }
}

event acceleration (i++; t <= time_pert)
{
  foreach_face(y)
    a.y[] += (t > time_pert? 0:
               -gpert*rhol*sq(dia)/sigma);
}

event logfile (i++; t <= t_end)
{
  double v1 = i == 0? v0: 0;
  foreach(reduction(+:v1))
    v1 += cr[] * pow(Delta, dimension);

  scalar pos[];
  position (f, pos, {0,1,0});
  double hmax = statsf(pos).max;
  position (cr, pos, {1,0,0});
  double rmax = statsf(pos).max;

  fprintf(stdout, "%d %g %g %g %g %g %g %d %g %g\n", 
      i, t, theta0, v0, v1, (v1 - v0)/(v0)*100, avgitr, g_count, hmax, rmax);
  fprintf(stderr, "%d %g %g %g %g %g %g %d %g %g\n", 
      i, t, theta0, v0, v1, (v1 - v0)/(v0)*100, avgitr, g_count, hmax, rmax);

  fflush(stdout);
}

event dump (t += 2.5)
{
  scalar * dumplist = {u,p,f,ibm,cr,ch};
  char name[80];
  sprintf(name, "dump-%g", t);
  dump (file = name, list = dumplist);
}

#if 1
event movie (t += 0.05; t <= t_end)
{
  char name[80];
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 30, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.299, tz = -2.748, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
  sprintf(name, "%g-f-%g.png", theta0, t);
  save (name);

  clear();
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "cr", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "cr", edges = true);
  sprintf(name, "%g-cr-%g.png", theta0, t);
  save (name);

  clear();
  cells (n = {0,1,0});
  //draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "ch", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "ch", edges = true);
  sprintf(name, "%g-ch-%g.png", theta0, t);
  save (name);
}
#endif

#if 1
event adapt (i++) {
  scalar ch1[], ibm1[];
  foreach() {
    ch1[] = f[]; 
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, ch1, u}, (double[]){1e-3, 1e-3, 3e-2, 3e-2, 3e-2}, 
    minlevel = min_level, maxlevel = max_level);
}
#endif

