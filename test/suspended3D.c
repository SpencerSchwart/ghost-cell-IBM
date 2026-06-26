/**
File for simulating a 3D suspended droplet on a flat, free-slip plate. The
droplet is perturbed in a way in which its first/second mode is excited.

This file uses Basilisk's built in contact.h to impose the boundary condition
and is used as a reference. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "vof.h"

#include "navier-stokes/perfs.h"
#include "view.h"

/** 
Simulation variables */

int maxlevel = 8;
int minlevel = 4;

double df = 1;          // droplet diameter (5 mm)
double theta0 = 90;  // contact angle
double v0 = 0;          // initial volume

double omegac = 23.664; // capillary frequency, sqrt(sigma/(rhol * D^3))
double gpert = 1.75;   // gravity perturbation to excite 1st mode (2 m/s^2)
double tpert = 0.56;    // time until perturbation ends (23.67 ms)

double tend = 40;
double tout = 0.1; 

FILE * outf;           // out file
double trestart = -1;  // time to restart (if given)

/**
Boundary conditions. */

vector h[];
h.t[bottom] = contact_angle (theta0*pi/180.);
h.r[bottom] = contact_angle (theta0*pi/180.);

u.n[bottom] = dirichlet(0);
u.t[bottom] = neumann(0);
u.r[bottom] = neumann(0);

u.n[top] = dirichlet(0);
p[top]   = dirichlet(0);
pf[top]  = dirichlet(0);


int main(int argc, char* argv[])
{
  size (4.);

  /**
  We shift the bottom boundary. */
  origin (-2, 0, -2);
  init_grid (1 << minlevel);
 
  /**
  Set the multigrid solver tolerances. */
  TOLERANCE    = 1e-6;

  /**
  We use fluid properties corresponding to water and air, normalized using the
  liquid density, surface tension, and the diameter of the droplet. */
  rho1 = 1;
  rho2 = 0.0012;
  mu1 = 1.69E-3; // Ohnesorge number
  mu2 = 3.04E-5;
  f.sigma = 1;

  if (argc > 1)
    trestart = atof(argv[1]);

  f.height = h;

  run();
}


event init (t = 0)
{
  char name[80];
  sprintf(name, "dumps/%g-dump-%g\n", theta0, trestart);
  if (!restore(name)) {

    /**
    Iteratively refine the initial mesh around the droplet. */
    scalar f1[];

    astats st;
    int count = 0;
    do {
      fraction (f1, - (sq(x) + sq(y) + sq(z) - sq(0.5*df)));
      st = adapt_wavelet ({f1}, (double[]){1e-5}, 
                          maxlevel = maxlevel, minlevel = minlevel);
      count++;
    } while ((st.nc || st.nf) && count < 20);

    fraction (f,  - (sq(x) + sq(y) + sq(z) - sq(0.5*df)));
 
    foreach(reduction(+:v0))
      v0 += f[]*cube(Delta);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "w");
    }
  }
  else {
    boundary(all);

    foreach(reduction(+:v0))
      v0 += f[]*cube(Delta);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "a");
    }
  }
}


/**
The perturbation is applied as an acceleration like gravity. Note the 
accleration is only applied to the liquid phase. */
event acceleration (i++)
{
  if (t <= tpert) {
    face vector av = a;
    foreach_face(y)
      av.y[] -= (f[] + f[0,-1])*0.5*gpert; 
  }
}


event logfile (i++; t <= tend)
{

  vector cl[];

  /**
  Calculate volume, maximum capillary number, and maximum droplet height. */

  double vreal = 0, ca = 0, hf = 0;
  foreach(reduction(+:vreal) reduction(max:ca) reduction(max:hf)) {

    foreach_dimension()
      cl.x[] = nodata;

    vreal += f[]*cube(Delta);

    double umag = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
    if (umag > ca)
      ca = umag;

    if (f[] > 0 && f[] < 1) {
      coord n = interface_normal(point, f), mp;
      double alpha = plane_alpha(f[], n);
      plane_area_center (n, alpha, &mp);

      mp.y = y + Delta*mp.y;

      if (mp.y > hf)
        hf = mp.y;
    }

  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  /**
  Get average radius of contact line along with some statistics like stddev 
  and centroid position. */

  foreach_boundary(bottom) {

    if (f[] > 1e-6 && f[] < 1 - 1e-6) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);

      coord v[12];
      int cp = facets (n, alpha, v, 1);

      coord pcl[2];
      int ccl = 0;
      for (int i = 0; i < cp; i++) {
        if (v[i].y <= -0.5 + 1e-10 && ccl < 2) {
          foreach_dimension()
            pcl[ccl].x = v[i].x;
          ccl++;
        }
      }

      if (ccl == 2) {
        coord cc = {x, y, z};
        foreach_dimension()
          cl.x[] = cc.x + Delta*0.5*(pcl[0].x + pcl[1].x);
      }
    }
  }

  double rcl = 0, clcx = 0, clcy = 0, clcz = 0;
  int count = 0;
  scalar scl[];
  foreach (reduction(+:rcl) reduction(+:clcx) reduction(+:clcy) reduction(+:clcz) reduction(+:count)) {
    scl[] = cl.x[] != nodata? sqrt(sq(cl.x[]) + sq(cl.y[]) + sq(cl.z[])): nodata;
    if (cl.x[] != nodata) {
      rcl  += scl[];
      clcx += cl.x[];
      clcy += cl.y[];
      clcz += cl.z[];
      count++;
    }
  }

  if (count > 0) {
    rcl  /= (double)count;
    clcx /= (double)count;
    clcy /= (double)count;
    clcz /= (double)count;
  }

  stats clstats = statsf(scl);

  if (pid() == 0 && outf) {
    fprintf(outf, "%d %g %g %g %g %g %g %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verror, rcl, ca, hf, clstats.stddev, clcx, clcy, clcz);
    fflush(outf);
  }
}


event dump (t += tout)
{
  scalar * dumplist = {u,p,f};
  char name[80];
  sprintf(name, "dumps/%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}


/**
Output images of the oscillating droplet and coordinates of the contact line. */
event movie (t += tout; t <= tend)
{
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 20, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.0, tz = -2.748, width = 700, height = 550);
  draw_vof ("ch", lw = 2, fc = {1, 0.62891, 0.62891});
  char name[80];
  sprintf(name, "imgs/%g-movie-%.5f.png", theta0, t);
  save (name);

  /**
  Output coordinates of the contact line. */
  sprintf(name, "data/%d-%g-contact-line-%g", pid(), theta0, t);
  FILE * fp = fopen(name, "w");

  vector cl[];
  foreach_boundary(bottom) {

    foreach_dimension()
      cl.x[] = nodata;

    if (f[] > 1e-6 && f[] < 1 - 1e-6) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);

      coord v[12];
      int cp = facets (n, alpha, v, 1);

      coord pcl[2];
      int ccl = 0;
      for (int i = 0; i < cp; i++) {
        if (v[i].y <= 1e-10 && ccl < 2) {
          foreach_dimension()
            pcl[ccl].x = v[i].x;
          ccl++;
        }
      }

      if (ccl == 2 && fp) {
        fprintf (fp, "%g %g %g\n%g %g %g\n\n", 
            x + pcl[0].x*Delta, y + pcl[0].y*Delta, z + pcl[0].z*Delta,
            x + pcl[1].x*Delta, y + pcl[1].y*Delta, z + pcl[1].z*Delta);
      }
    }
  }

  if (fp)
    fclose(fp);
}


/**
We trick the adapt function to keep the liquid interface at the maximum level
of refinement. */
event adapt (i++) 
{
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1, u}, (double[]){1e-3, 5e-2, 5e-2, 5e-2}, 
    minlevel = minlevel, maxlevel = maxlevel);
}

