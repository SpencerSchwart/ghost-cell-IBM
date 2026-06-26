/**
File for simulating a 3D suspended droplet on a flat, free-slip plate. The
droplet is perturbed in a way in which its first/second mode is excited. */

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
Boundary conditions */

// Free-slip wall
u.n[immersed] = dirichlet(0);
u.t[immersed] = neumann(0);
u.r[immersed] = neumann(0);

u.n[top] = neumann(0);
p  [top] = dirichlet(0);
pf [top] = dirichlet(0);

u.n[right] = neumann(0);
u.n[front] = neumann(0);

int main(int argc, char* argv[])
{
  size (4.);

  /**
  We shift the bottom boundary. */
  origin (-2, -0.254, -2);
  init_grid (1 << minlevel);
 
  /**
  Set the multigrid solver tolerances. */
  TOLERANCE = 1e-6;

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

  f.wetting.theta_s = theta0;

  run();
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
  sprintf(name, "dumps/%g-dump-%g\n", theta0, trestart);
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

    fraction (f,  - (sq(x) + sq(y) + sq(z) - sq(0.5*df)));
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

  /**
  Calculate volume, maximum capillary number, and maximum droplet height. */

  double vreal = 0, ca = 0, hf = 0;
  foreach(reduction(+:vreal) reduction(max:ca) reduction(max:hf)) {

    vreal += f[]*cube(Delta);

    if (gc[] && cs[] == 1) {
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
  Get average radius of contact line along with some statistics like stddev 
  and centroid position. */
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
        i, t, theta0, v0, vreal, verror, rcl, ca, hf, clstats.stddev, clcx, clcy, clcz);
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
Output images of the oscillating droplet and coordinates of the contact line. */
event movie (t += tout; t <= tend)
{
  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 20, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.0, tz = -2.748, width = 700, height = 550);
  draw_vof("cs", "fs", fc = {0.796875, 0.796875, 0.796875});
  draw_vof("cs", "fs", edges = true);
  draw_vof ("ch", lw = 2, fc = {1, 0.62891, 0.62891});
  char name[80];
  sprintf(name, "imgs/%g-movie-%.5f.png", theta0, t);
  save (name);

  sprintf(name, "data/%d-%g-contact-line-%g", pid(), theta0, t);
  FILE * fp = fopen(name, "w");
  output_contact_line (ch, cs, fs, fp);
  fclose(fp);
}


/**
We trick the adapt function to keep the liquid interface at the maximum level
of refinement. Note we do not do that with the solid interface. */
event adapt (i++)
{
  scalar f1[];
  foreach()
    f1[] = f[];

  adapt_wavelet ({f1, u}, (double[]){1e-3, 5e-2, 5e-2, 5e-2}, 
    maxlevel = maxlevel, minlevel = minlevel);

  solid_domain (cs, fs);
}

