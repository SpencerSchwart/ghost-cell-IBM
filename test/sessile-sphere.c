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
double Oh = 0.1;       // Ohnesorge number
double v0 = 0;         // initial volume

double tend = 4;
double tout = 0.25;

FILE * outf;           // out file
double trestart = -1;  // time to restart (if given)

/**
The droplet is initalized so that it corresponds to the equilibirum shape of 
contact angle = 90 degrees, meaning the droplet will spread when CA < 90 and
retract when CA > 90. 

TODO: it is now initalized a bit higher to avoid droplet launch off at CA = 165. */

coord xd0 = {0, sqrt(2.)/1.5, 0};

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
  size (4.);

  origin (0, -2, 0);

  init_grid(1 << (minlevel));
  
  TOLERANCE    = 1e-4;
  TOLERANCE_MU = 1e-3;

  /**
  We use a constant fluid properties. */
  rho1 = rho2 = 1;
  f.sigma = 0.1;
  mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

  if (argc > 1)
    theta0 = atof(argv[1]);
  if (argc > 2)
    tend = atof(argv[2]);
  if (argc > 3)
    maxlevel = atoi(argv[3]);
  if (argc > 4)
    trestart = atof(argv[4]);

  /**
  We vary the contact_angle if theta0 != 0. */
  if (theta0) {
    f.wetting.theta_s = theta0;
    run();
  } 
  else {
     double angles[7] = {15,30,60,90.05,120,150,165};
     for (int i = 0; i < 7; i++) {
       theta0 = angles[i];
       tend = theta0 == 15? 35: theta0 == 30? 30: 15;
       f.wetting.theta_s = theta0;
       run();
     }
  }
}


event init (t = 0)
{
  char name[80];
  sprintf(name, "dumps/%g-dump-%g\n", theta0, trestart);
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
  }
  else {
    solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
    boundary(all);
  }
}


event logfile (i++; t <= tend)
{
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). We check this every 10 iterations to
  reduce computational costs associated with calculating the curvature. */
 
  if (i % 10 == 0) {
     scalar kappa[];
     curvature (ch, kappa);
     foreach()
       if (cs[] < 1.)
         kappa[] = nodata;
     if (i > 10 && statsf (kappa).stddev < 1e-3)
       return true;
  }

  /**
  Calculate volume, and maximum capillary number. */

  double vreal = 0, ca = 0;
  foreach(reduction(+:vreal) reduction(max:ca)) {

    vreal += f[]*cube(Delta);

    if (gc[] && cs[] == 1) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
      if (umag > ca)
        ca = umag;
    }
  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  /**
  Calculate maximum droplet height and wetted contact line radius.
  To account for the curvature of the solid substrate, we calculate the angular
  position of the contact line as well. 

  TODO: adapt contact line position to 3D better. */

  double hf = 0, anglecl = -HUGE;
  foreach_boundary(back, reduction(max:hf) reduction(max:anglecl)) {
    if (cs[] == 1 && f[] > 0 && f[] < 1) {
      coord n = interface_normal(point, f), mp;
      double alpha = plane_alpha(f[], n);
      plane_area_center (n, alpha, &mp);
      mp = local_to_global_coord(point, mp);
      if (mp.y > hf)
        hf = mp.y;
    }

    if (extra[]) {
      coord n = interface_normal(point, ch), mp;
      double alpha = plane_alpha(ch[], n);

      plane_area_center (n, alpha, &mp);
      mp = local_to_global_coord(point, mp);

      double angle = atan2(mp.y, mp.x);
      if (angle > anglecl)
        anglecl = angle;
    }
  }

  double rcl = anglecl*(ds*0.5);

  /**
  If we are restarting the simulation, append to the existing out file. */
  if (i == 0) {
    char name[80];
    sprintf(name, "outs/out-%g", theta0);
    outf = trestart > -1? fopen(name, "a"): fopen(name, "w");
  }

  fprintf(outf, "%d %g %g %g %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verror, rcl, ca, hf, anglecl*180/pi);
  fflush(outf);
}


event dump (t += tout)
{
  scalar * dumplist = {u,p,f,ch,cs,contact_angle,extra,gginter};
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
  sprintf(name, "imgs/%g-chleft-%g.png", theta0, t);
  save (name);

  view (camera="right", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells ();
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  squares ("ch");
  sprintf(name, "imgs/%g-chright-%g.png", theta0, t);
  save (name);

  view (camera="back", fov = fov, bg = {1,1,1}, width = 700, height = 550, ty = ty);
  cells (n = {0,1,0});
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  squares ("ch");
  sprintf(name, "imgs/%g-chback-%g.png", theta0, t);
  save (name);
}


event end (t = end)
{
  /** 
  At equilibrium we output the (almost constant) radius, volume, maximum velocity and time.*/
  char name[80];
  sprintf(name, "results/results-%g", theta0);
  FILE * fp = fopen (name, "w");

  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1.)
        kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 4*statsf(f).sum;

  fprintf (fp, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
     N, theta0, R, s.stddev, V, normf(u.x).max, t);
  fflush(fp);
  fclose(fp);

  scalar * dumplist = {u,p,f,ch,cs,contact_angle,extra,gginter};
  sprintf(name, "dumps/%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}


event adapt (i++) 
{
  scalar f1[];
  foreach() 
    f1[] = ch[];

  adapt_wavelet ({cs, f1}, (double[]){1e-3, 1e-3}, 
    minlevel = minlevel, maxlevel = maxlevel);

  solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
  fractions_cleanup(cs, fs);
}

