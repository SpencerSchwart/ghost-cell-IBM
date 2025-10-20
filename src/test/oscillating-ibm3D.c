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
#define T 2
#define R0 0.2

double Oh = 0.05; // Ohnesorge number
double theta0=150, volume_vof_init;

#define MAXLEVEL 7
#define MINLEVEL 3

u.t[immersed] = dirichlet(0);
u.n[immersed] = neumann(0);
u.r[immersed] = neumann(0);

u.n[left] = dirichlet(0);
u.t[left] = neumann(0);
u.n[back] = dirichlet(0);
u.t[back] = neumann(0);

f[left] = neumann(0);
f[back] = neumann(0);

p[left] = neumann(0);
p[back] = neumann(0);

pf[left] = neumann(0);
pf[back] = neumann(0);

p[right]  = dirichlet(0);
pf[right] = dirichlet(0);

int main(int argc, char * argv[])
{
  theta0 = 70;

  if (argc > 1)
    theta0 = atof (argv[1]);

  size (L0);
  init_grid(1 << MINLEVEL);
  origin (0, -0.26, 0);

  //DT = 0.01;

  TOLERANCE = 1e-5;
#if 0
  rho1 = rho2 = 1;
  f.sigma = 0.1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * 2.*R0);
#else
  mu1 = 0.1; mu2 = 0.1;
  rho1= 0.1; rho2 = 0.1;
  f.sigma = 1.;
#endif
#if 1
  for (theta0 = 30; theta0 <= 150; theta0 += 30) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
#else
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
#endif
}

void solid_domain (scalar s, face vector sf)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y + 0.00725;
    //phi[] = y + 0.00725;
  boundary ({phi});
  fractions (phi, s, sf);
  boundary({s, sf});
}

double v0;
event init (t = 0)
{
  if(!restore("70-dump-1.55")) {  
    
    scalar ibmFake[], fFake[];
    face vector ibmfFake[];

    astats st; int count = 0;
    do {
      solid_domain(ibmFake, ibmfFake);
      fraction (fFake, - (sq(x) + sq(y) + sq(z) - sq(R0)));
      st = adapt_wavelet ({ibmFake, fFake}, (double[]){1e-5, 1e-5}, MAXLEVEL, MINLEVEL);
      count++;
    } while ((st.nc || st.nf) && count < 30);

    solid_domain (ibm, ibmf);
    fraction (f, - (sq(x) + sq(y) + sq(z) - sq(R0)));
    fraction (ch, - (sq(x) + sq(y) + sq(z) - sq(R0)));
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

event logfile (i++; t <= T)
{
  double v1 = i == 0? v0: 0;
  foreach(reduction(+:v1))
    v1 += cr[] * pow(Delta, dimension);

  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4.*statsf(cr).sum;

  if (statsf (kappa).stddev < 1e-3)
    return true;

  scalar pos[];
  position (f, pos, {0,1,0});
  double hmax = statsf(pos).max;
  position (cr, pos, {1,0,0});
  double rmax = statsf(pos).max;

  fprintf(stdout, "%d %g %g %g %g %g %g %d %g %g %g %g\n", 
      i, t, theta0, v0, v1, (v1 - v0)/(v0)*100, avgitr, g_count, hmax, rmax, R, V);
  fprintf(stderr, "%d %g %g %g %g %g %g %d %g %g %g %g\n", 
      i, t, theta0, v0, v1, (v1 - v0)/(v0)*100, avgitr, g_count, hmax, rmax, R, V);

  fflush(stdout);
}

#if DUMP
event dump (t += 1)
{
  scalar * dumplist = {u,p,f,ibm,cr,ch,gf0};
  char name[80];
  sprintf(name, "%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}
#endif

#if 1
event movie (t += 0.05; t <= T)
{
  char name[80];

#if 0
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 30, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.299, tz = -2.748, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
  mirror (n = {1,0,0}) {
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
    }
  mirror (n = {0,0,1}) {
  cells (n = {0,1,0});
  draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
    mirror (n = {1,0,0}) {
    cells (n = {0,1,0});
    draw_vof (c = "ibm", s = "ibmf", filled = -1, fc = {0.8,0.8,0.8});
    draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
    draw_vof (c= "f", edges = true);
    }
  }
  sprintf(name, "%g-movie-%g.png", theta0, t);
  save (name);
#endif

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

event end (t = end)
{
  FILE * fp = fopen ("results", "a");
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
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
  adapt_wavelet ({f1}, (double[]){1e-3}, 
    maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

