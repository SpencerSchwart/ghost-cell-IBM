/**
Droplet impacting on a 3D sphere. */

#include "grid/octree.h"
#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

#include "view.h"
#include "tag.h"
#include "navier-stokes/perfs.h"

/**
Simulation variables. */

int maxlevel = 12;
int minlevel = 6;

double xd0 = 1.9;       // initial drop position
double ud0 = 1;         // initial drop velocity

double df = 1;          // drop diameter
double ds = 2.7;        // sphere diameter

double theta0 = 120;    // contact angle

double gravity = 1;
double v0 = 0;          // initial volume

double tout = 0.05;
double tend = 12;
double trestart = 1000;
double uemax = 1e-3;
double cflg = 0.5;

FILE * outf;

/**
Boundary conditions. */

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

u.n[right] = neumann(0);
p  [right] = dirichlet(0);
pf [right] = dirichlet(0);

u.n[left] = neumann(0);
p  [left] = dirichlet(0);
pf [left] = dirichlet(0);


int main(int argc, char* argv[])
{
  double l0 = 50;

  /**
  Read the input parameters. */
  if (argc > 1)
    l0 = atof(argv[1]);
  if (argc > 2)
    maxlevel = atoi(argv[2]);
  if (argc > 3)
    minlevel = atoi(argv[3]);
  if (argc > 4)
    rho1 = atof(argv[4]);
  if (argc > 5)
    rho2 = atof(argv[5]);
  if (argc > 6)
    mu1 = atof(argv[6]);
  if (argc > 7)
    mu2 = atof(argv[7]);
  if (argc > 8)
    f.sigma = atof(argv[8]);
  if (argc > 9)
    gravity = atof(argv[9]);
  if (argc > 10)
    ud0 = atof(argv[10]);
  if (argc > 11)
    f.wetting.theta_s = atof(argv[11]);
  if (argc > 12)
    xd0 = atof(argv[12]);
  if (argc > 13)
    df = atof(argv[13]);
  if (argc > 14)
    ds = atof(argv[14]);
  if (argc > 15)
    tend = atof(argv[15]);
  if (argc > 16)
    tout = atof(argv[16]);
  if (argc > 17)
    trestart = atof(argv[17]);
  if (argc > 18)
    DT = atof(argv[18]);
  if (argc > 19)
    cflg = atof(argv[19]);

  /**
  Set multigrid tolerances. */
  TOLERANCE    = 1e-4;

  double hmin = l0/(1 << maxlevel);

  size(l0);
  origin (-l0/2. -0.1*hmin, 0);
  init_grid (1 << minlevel);

  run();
}

event init (t = 0)
{
  CFL = cflg;

  char name[80];
  sprintf(name, "dumps/dump-%g", trestart);
  if (!restore(name)) {
    scalar cs1[], f1[];
    face vector fs1[];

    /**
    Iteratively refine the initial mesh around the sphere and droplet. */
    int count = 0;
    astats st;
    do {
      solid (cs1, fs1, (sq(x) + sq(y) + sq(z) - sq(ds*0.5)));
      fraction (f1, - (sq(x - xd0) + sq(y) + sq(z) - sq(df*0.5)));
      st = adapt_wavelet ({cs1, f1}, (double[]){1e-5, 1e-5}, maxlevel, minlevel);
      count++;
    } while ((st.nc || st.nf) && count < 40);

    solid (cs, fs, (sq(x) + sq(y) + sq(z) - sq(ds*0.5)));
    fractions_cleanup(cs, fs, 1e-3);

    fraction(f,  - (sq(x - xd0) + sq(y) + sq(z) - sq(df*0.5)));
    fraction(ch, - (sq(x - xd0) + sq(y) + sq(z) - sq(df*0.5)));

    foreach()
      u.x[] = -ud0*f[]*cs[];

    v0 = real_volume(f, cs, fs);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "w");
    }
  }
  else {
    solid (cs, fs, (sq(x) + sq(y) + sq(z) - sq(ds*0.5)));
    boundary({cs, fs});
    fractions_cleanup(cs, fs, 1e-3);
    restriction({cs, fs});

    boundary(all);

    v0 = real_volume(f, cs, fs);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "a");
    }
  }
}


/**
Apply gravity to the accleration field, where left is down. */
event acceleration (i++) 
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= gravity;
}


event logfile (i++; t <= tend)
{
  /**
  Calculate the volume loss/gain error. */
  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += f[]*cube(Delta);

  double verr = i == 0 && v0? 0: (vreal - v0)/v0 * 100;

  /**
  Calculate the droplet height at the impact location. */
  double hf = 0;
  foreach_boundary(bottom, reduction(max:hf)) {
    if (f[] > 1e-6 && f[] < 1- 1e-6 && cs[] == 1) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double htemp = x + p.x*Delta - 0.5*ds;
      if (htemp > hf)
        hf = htemp;
    }
  }

  /**
  Get average radius and stddev of the contact line. */
  vector cl[];
  measure_contact_line(cl, ch, cs, fs);

  scalar scl[];
  foreach() {
    if (cl.x[] == nodata)
      scl[] = nodata;
    else {
      coord a;
      foreach_dimension()
        a.x = cl.x[];

      /**
      Calculate the angle between the postion vector a and the central axis of
      the sphere. */

      scl[] = dot_product_angle(a, (coord){1,0,0});
    }
  }
  
  stats clstats = statsf(scl);
  norm clnorm = normf(scl);

  if (pid() == 0) {
    fprintf(outf, "%d %g %g %g %g %g %g %g %g %g\n", 
      i, t, theta0, v0, vreal, verr, hf, clnorm.avg*ds, clnorm.avg*180./pi, clstats.stddev);
    fflush(outf);
  }
}


event dump (t += tout)
{
  scalar * dumplist = {u,p,f,cs,ch,extra,contact_angle,gginter};
  char name[80];
  sprintf(name, "dumps/dump-%g", t);
  dump (file = name, list = dumplist);
}


event movie (t += tout) 
{
  /**
  We output images from a few different views. */
  char name[80];
  sprintf(name, "imgs/%g-back-%.5f.png", theta0, t);
  view(camera="back", fov = 6, tx = -0.025, width = 1000, height = 1000, bg={1,1,1});
  clear();
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  save(name);

  sprintf(name, "imgs/%g-crfront-%.5f.png", theta0, t);
  view(camera="front", fov = 6, tx = -0.025, width = 1000, height = 1000, bg={1,1,1});
  clear();
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  save(name);

#if 0
  sprintf(name, "imgs/%g-crfrontu-%.5f.png", theta0, t);
  view(camera="front", fov = 6, tx = -0.025, width = 1000, height = 1000, bg={1,1,1});
  clear();
  draw_vof ("ch", lw = 2, color="u.x", min=-1, max=1, linear=true);
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  mirror (n = {0,-1,0}) {
    draw_vof ("ch", lw = 2, color="u.x", min=-1, max=1, linear=true);
    draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  }
  save(name);

  sprintf(name, "imgs/%g-crisou-%.5f.png", theta0, t);
  view(camera="iso", fov = 6, tx = -0.025, width = 1000, height = 1000, bg = {1,1,1});
  clear();
  draw_vof ("ch", lw = 2, color="u.x", min=-1, max=1, linear=true);
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  mirror (n = {0,0,-1}) {
    draw_vof ("ch", lw = 2, color="u.x", min=-1, max=1, linear=true);
    draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  }
  save(name);
#endif

  /**
  We output the contact line coordinates. */
  sprintf(name, "data/%d-%g-contact-line-%g", pid(), theta0, t);
  FILE * fp = fopen(name, "w");
  output_contact_line (ch, cs, fs, fp);
  fclose(fp);

  /**
  Lastly, we output the VOF reconstructed surfaces. */
  sprintf(name, "data/%d-%g-vof-%g", pid(), theta0, t);
  fp = fopen(name, "w");
  output_facets_contact (f, ch, cs, fp);
  fclose(fp);

  sprintf(name, "data/%d-%g-cs-%g", pid(), theta0, t);
  fp = fopen(name, "w");
  output_facets (cs, fp, s = fs);
  fclose(fp);
}


/**
The mesh is adapted only every four iterations to increase speed. */
event adapt (i++)
{
  if (i % 4 == 0) {
    scalar f1[];
    foreach() 
      f1[] = f[];
    
    adapt_wavelet ({cs, f1, u}, (double[]){5e-1, 1e-3, 1e-2, 1e-2, 1e-2}, 
      maxlevel = maxlevel, minlevel = minlevel);

    solid (cs, fs, (sq(x) + sq(y) + sq(z) - sq(ds*0.5)));
    fractions_cleanup(cs, fs, 1e-3);
  }
}

