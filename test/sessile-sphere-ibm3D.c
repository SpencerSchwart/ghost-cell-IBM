/**
A 3D sessile droplet spreading on a sphereical particle until it reaches an
equilibirum position, which is based on the contact angle. */

#include "grid/octree.h"
#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

#include "navier-stokes/perfs.h"
#include "view.h"

/** 
Simulation variables */

int maxlevel = 7;
int minlevel = 3;

double df = 1;         // droplet diameter
double ds = 1;         // sphere diameter
double theta0 = 0;     // contact angle
double Oh = 0.04;       // Ohnesorge number
double v0 = 0;         // initial volume

double tend = 4;
double tout = 0.25;

FILE * outf;           // out file
double trestart = -1;  // time to restart (if given)

/**
The droplet is initalized so that it corresponds to the equilibirum shape of 
contact angle = 90 degrees, meaning the droplet will spread when CA < 90 and
retract when CA > 90. */

coord xd0 = {0, 0.707106781, 0}; // sqrt(2)/2

/**
Boundary conditions. */

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

u.n[bottom] = neumann(0);
p  [bottom] = dirichlet(0);
pf [bottom] = dirichlet(0);

int main(int argc, char* argv[])
{
  if (argc > 1)
    theta0 = atof(argv[1]);
  if (argc > 2)
    tend = atof(argv[2]);
  if (argc > 3)
    maxlevel = atoi(argv[3]);
  if (argc > 4)
    Oh = atof(argv[4]);
  if (argc > 5)
    trestart = atof(argv[5]);

  size (4.);

  origin(-2, -2, -2);
  init_grid (1 << minlevel);
  
  TOLERANCE = 1e-3;

  /**
  We use constant fluid properties. */
  rho1 = rho2 = 1;
  f.sigma = 1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

  /**
  We vary the contact_angle if theta0 != 0. */
  if (theta0) {
    f.wetting.theta_s = theta0;
    run();
  } 
  else {
     double angles[7] = {15,30,60,90,120,150,165};
     for (int i = 0; i < 7; i++) {
       theta0 = angles[i];
       tend = theta0 == 15? 14: theta0 == 30? 10: 10;

       if (theta0 >= 150)
         Oh = 0.16;

       mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

       f.wetting.theta_s = theta0;
       run();
     }
  }
}


event init (t = 0)
{
  char name[80];
  sprintf(name, "dumps/%g-dump-%g", theta0, trestart);
  if (!restore(name)) {

    /**
    Iteratively refine the initial mesh around the immersed boundary and droplet.
    Here, the immersed boundary / solid is a sphere centered at the origin, which
    has the same diamter as the droplet. */

    scalar cs1[], f1[];
    face vector fs1[];

    astats st;
    int count = 0;
    do {
      solid (cs1, fs1, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
      fraction (f1, - (sq(x - xd0.x) + sq(y - xd0.y) + sq(z - xd0.z) - sq(0.5*df)));
      st = adapt_wavelet ({cs1, f1}, (double[]){1e-5, 1e-5}, 
                          maxlevel = maxlevel, minlevel = minlevel);
      count++;
    } while ((st.nc || st.nf) && count < 20);

    /**
    Initalize the geometry and both volume fraction fields, f and ch. Also
    calculate the initial volue v0, which is the portion of f outside cs and fs. */

    solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
    event("update_metric");

    fraction (f,  - (sq(x - xd0.x) + sq(y - xd0.y) + sq(z - xd0.z) - sq(0.5*df)));
    fraction (ch, - (sq(x - xd0.x) + sq(y - xd0.y) + sq(z - xd0.z) - sq(0.5*df)));
  
    v0 = real_volume(f, cs, fs);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "w");
    }
  }
  else {
    solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
    boundary(all);

    foreach(reduction(+:v0))
      v0 += f[]*cube(Delta);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "a");
    }
  }
}


event logfile (i++; t <= tend)
{
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). We check this every 10 iterations to
  reduce computational costs associated with calculating the curvature. */
 
  if (i % 5 == 0) {
     scalar kappa[];
     curvature (ch, kappa);
     foreach()
       if (cs[] < 1. || f[] > cs[] - 1e-6)
         kappa[] = nodata;
     if (i > 10 && statsf (kappa).stddev < 1e-3)
       return true;
  }

  /**
  Calculate volume, maximum capillary number, and maximum droplet height. */

  double vreal = 0, ca = 0, hf = 0;
  foreach(reduction(+:vreal) reduction(max:ca) reduction(max:hf)) {

    vreal += f[]*cube(Delta);

    if (gc[] && cs[] == 1 && f[]) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
      if (umag > ca)
        ca = umag;
    }

    if (cs[] == 1 && f[] > 0 && f[] < 1) {
      coord n = interface_normal(point, f), mp;
      double alpha = plane_alpha(f[], n);
      plane_area_center (n, alpha, &mp);
      mp = local_to_global_coord(point, mp);

      if (mp.y > hf)
        hf = mp.y;
    }
  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  /**
  Get average radius and stddev of the contact line. */
  vector cl[];
  double rcl = measure_contact_line(cl, ch, cs, fs);

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
  stats2 extras =  statsf2(extra, val = 0);

  if (pid() == 0) {
    fprintf(outf, "%d %g %g %g %g %g %g %g %g %g %g %g %g\n", 
        i, t, theta0, v0, vreal, verror, rcl, clstats.stddev, ca, hf, 
        extras.sum, clnorm.avg*ds, clnorm.avg*180./pi);
    fflush(outf);
  }
}


event dump (t += tout)
{
  scalar * dumplist = {u,p,f,ch,cs,contact_angle,extra};
  char name[80];
  sprintf(name, "dumps/%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}


event movie (t += tout; t <= tend)
{
  char name[80];
  double fov = 18, ty = -0.2;

  view (camera="left", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  sprintf(name, "imgs/%g-chleft-%.5f.png", theta0, t);
  save (name);

  view (camera="right", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells ();
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  squares ("ch");
  sprintf(name, "imgs/%g-chright-%.5f.png", theta0, t);
  save (name);

  view (camera="back", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  squares ("ch");
  sprintf(name, "imgs/%g-chback-%.5f.png", theta0, t);
  save (name);
}


event end (t = end)
{
  if (pid() == 0)
    fclose(outf);

  /** 
  At equilibrium we output the (almost constant) radius, volume, maximum velocity and time.*/
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1. || f[] > cs[] - 1e-6)
        kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;

  char name[80];
  sprintf(name, "results/results-%g", theta0);
  FILE * fp = fopen (name, "w");
  fprintf (fp, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
     N, theta0, R, s.stddev, V, normf(u.x).max, t);

  fclose (fp);

  /**
  We output the dumpfile of the last timestep. */
  scalar * dumplist = {u,p,f,ch,cs,contact_angle,extra};
  sprintf(name, "dumps/%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);

  /**
  We output the contact line coordinates. */
  sprintf(name, "data/%d-%g-contact-line-%g", pid(), theta0, t);
  fp = fopen(name, "w");
  output_contact_line (ch, cs, fs, fp);
  fclose(fp);

  /**
  Lastly, we output the VOF reconstructed surface. */
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
We trick the adaptation function to always keep the liquid interface fully
refined. */

event adapt (i++) 
{
  scalar f1[];
  foreach()
    f1[] = f[];
  
  adapt_wavelet ({cs, f1}, (double[]){1e-3, 1e-3}, 
    maxlevel = maxlevel, minlevel = minlevel);

  solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
  fractions_cleanup(cs, fs, 1e-3);
}

