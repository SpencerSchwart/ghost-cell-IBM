/**
A 3D sessile droplet oscillating on a flat, horizontal plate until it reaches an
equilibirum position, which is based on the contact angle. The oscillation occurs
from low Ohnesorge number and initalizing the droplet at 90 degrees instead of
its equilibirum contact angle (60 or 120). */

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

int maxlevel = 8;
int minlevel = 4;

double df = 1;         // droplet diameter
double theta0 = 0;     // contact angle
double Oh = 0.004;   // Ohnesorge number
double v0 = 0;         // initial volume

double tend = 20;
double tout = 0.25;    // 

FILE * outf;           // out file
double trestart = -1;  // time to restart (if given)

/** 
Boundary conditions */

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

u.n[top] = neumann(0);
p  [top] = dirichlet(0);
pf [top] = dirichlet(0);

int main(int argc, char* argv[])
{
  size (5.);

  /**
  We shift the bottom boundary. */
  origin (-2.5, -0.254, -2.5);
  init_grid (1 << minlevel);
  
  TOLERANCE = 1e-6;

  /**
  We use a constant fluid properties. */
  rho1 = rho2 = 1;
  f.sigma = 1;
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
     double angles[2] = {110,70};
     for (int i = 0; i < 2; i++) {
       theta0 = angles[i];
       f.wetting.theta_s = theta0;
       run();
     }
  }
}

/**
The geoemtry is a flat, horizontal plate at y = 0. */
void solid_domain (scalar cs, face vector fs)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y;
  boundary ({phi});
  fractions (phi, cs, fs);
}

event init (t = 0)
{
  char name[80];
  sprintf(name, "dumps/%g-dump-%g", theta0, trestart);
  if (!restore(name)) {

    /**
    Iteratively refine the initial mesh around the plate and droplet. */
    scalar cs1[], f1[];
    face vector fs1[];

    astats st;
    int count = 0;
    do {
      solid_domain (cs1, fs1);
      fraction (f1, - (sq(x) + sq(y) + sq(z) - sq(0.5*df)));
      st = adapt_wavelet ({cs1, f1}, (double[]){1e-5,1e-5}, 
                          maxlevel = maxlevel, minlevel = minlevel);
      count++;
    } while ((st.nc || st.nf) && count < 20);

    /**
    Initalize the geometry and both volume fraction fields, f and ch. Also
    calculate the initial volue v0, which is the portion of f outside cs and fs. */
    solid_domain (cs, fs);
    event("update_metric");
    fraction (f, - (sq(x) + sq(y) + sq(z) - sq(0.5*df)));
    fraction (ch, - (sq(x) + sq(y) + sq(z) - sq(0.5*df)));
  
    v0 = real_volume(f, cs, fs);

    if (pid() == 0) {
      sprintf(name, "outs/out-%g", theta0);
      outf = fopen(name, "w");
    }
  }
  else {
    solid_domain (cs, fs);
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
  Calculate volume, maximum drolet height, and maximum capillary number. */

  double vreal = 0, hf = 0, ca = 0;
  foreach(reduction(+:vreal) reduction(max:hf) reduction(max:ca)) {

    vreal += f[]*cube(Delta);

    if (cs[] == 1 && f[] > 0 && f[] < 1) {
      coord n = interface_normal(point, f), mp;
      double alpha = plane_alpha(f[], n);
      plane_area_center (n, alpha, &mp);
      mp = local_to_global_coord(point, mp);
      if (mp.y > hf)
        hf = mp.y;
    }

    if (gc[] && cs[] == 1 && f[]) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
      if (umag > ca)
        ca = umag;
    }
  }

  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  /**
  Get average radius and stddev of contact line. Also calculate the centroid of
  the contact line. */
  vector cl[];
  double rcl = measure_contact_line(cl, ch, cs, fs);

  double clcx = 0, clcy = 0, clcz = 0;
  int count = 0;
  scalar scl[];
  foreach(reduction(+:clcx) reduction(+:clcy) reduction(+:clcz) reduction(+:count)) {
    scl[] = cl.x[] != nodata? sqrt(sq(cl.x[]) + sq(cl.y[]) + sq(cl.z[])): nodata;
    if (cl.x[] != nodata) {
      clcx += cl.x[];
      clcy += cl.y[];
      clcz += cl.z[];
      count++;
    }
  }

  if (count > 0) {
    clcx /= (double)count;
    clcy /= (double)count;
    clcz /= (double)count;
  }

  stats clstats = statsf(scl);

  if (pid() == 0 && outf) {
    fprintf(outf, "%d %g %g %g %g %g %g %g %g %g %g %g %g\n", 
        i, t, theta0, v0, vreal, verror, rcl, clstats.stddev, ca, hf, clcx, clcy, clcz);
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


/**
Output images of the spreading/receding droplet. */
event movie (t += tout; t <= tend)
{
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 20, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.0, tz = -2.748, width = 700, height = 550);
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  draw_vof("cs", "fs", edges = true);
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});

  char name[80];
  sprintf(name, "imgs/%g-movie-%.5f.png", theta0, t);
  save (name);

  /**
  We output the contact line coordinates. */
  sprintf(name, "data/%d-%g-contact-line-%g", pid(), theta0, t);
  FILE * fp = fopen(name, "w");
  output_contact_line (ch, cs, fs, fp);
  fclose(fp);
}


event end (t = end)
{
  if (pid() == 0)
    fclose(outf);

  /** 
  At equilibrium we output the (almost constant) radius, volume, maximum velocity and time.*/
  char name[80];
  sprintf(name, "results/results-%g", theta0);
  FILE * fp = fopen (name, "w");

  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1. || f[] > cs[] - 1e-6)
        kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (fp, "%d %g %.5g %.3g %.5g %.3g %.5g\n",
     N, theta0, R, s.stddev, V, normf(u.x).max, t);

  fclose (fp);

  /**
  We output the dumpfile of the last timestep. */
  scalar * dumplist = {u,p,f,ch,cs,contact_angle};
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

event adapt (i++) 
{
  scalar f1[];
  foreach()
    f1[] = f[];

  adapt_wavelet ({cs, f1, u}, (double[]){1e-3, 1e-3, 5e-2, 5e-2, 5e-2}, 
    maxlevel = maxlevel, minlevel = minlevel);

  solid_domain (cs, fs);
}

