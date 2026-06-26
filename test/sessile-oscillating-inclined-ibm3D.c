/**
A 3D sessile droplet spreading on a flat, inclined plate until it reaches an
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

int maxlevel = 8;
int minlevel = 4;

double df = 1;         // droplet diameter
double theta0 = 0;     // contact angle
double Oh = 0.004;       // Ohnesorge number
double v0 = 0;         // initial volume

double tend = 40;
double tout = 0.25;      

FILE * outf;           // out file
double trestart = -1;  // time to restart (if given)

/**
Boundary conditions. */

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
  The origin if shifted down a bit (-0.05) so we can simply set the solid plate to be y = x. */
  origin (-2.5,-2.50001, -2.5);

  init_grid(1 << minlevel);
  
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
     double angles[2] = {110, 70};
     for (int i = 0; i < 2; i++) {
       theta0 = angles[i];
       f.wetting.theta_s = theta0;
       run();
     }
  }
}


/**
The geoemtry is a flat, inclined plate. */
void solid_domain (scalar c, face vector cf)
{
  vertex scalar phi[];
  foreach_vertex()
  phi[] = y + x + z;
  boundary ({phi});
  fractions (phi, c, cf);
  boundary({c, cf});
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

      double hfr = dot_product(mp, (coord){1./sqrt(3),1./sqrt(3),1./sqrt(3)});

      if (hfr > hf)
        hf = hfr;
    }

  }
  double verror = i == 0? 0: (vreal - v0)/(v0 + SEPS)* 100;
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
  scalar * dumplist = {u,p,f,ch,cs,contact_angle,extra,gginter};
  char name[80];
  sprintf(name, "dumps/%g-dump-%g", theta0, t);
  dump (file = name, list = dumplist);
}


/**
Output images of the spreading/receding droplet. */

event movie (t += tout; t <= tend)
{
  char name[80];

  view (camera="iso", fov = 30, near = 0.01, far = 1000, bg = {1,1,1},
        tx = 0, ty = 0.299, tz = -2.748, width = 1000, height = 1000);
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", edges = true);
  sprintf(name, "imgs/%g-chiso-%.5f.png", theta0, t);
  save (name);

  float q1[4], q2[4], q3[4];
  gl_axis_to_quat ((float[]){1,0,0}, -pi/10., q1);
  gl_axis_to_quat ((float[]){0,1,0}, -pi/4., q2);
  gl_add_quats(q1, q2, q3);

  view (quat=q3, fov = 20, bg = {1,1,1}, width = 1000, height = 1000);
  clear();
  draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
  draw_vof("cs", "fs", fc={0.796875, 0.796875, 0.796875});
  draw_vof("cs", "fs", edges = true);
  sprintf(name, "imgs/%g-chiso2-%.5f.png", theta0, t);
  save (name);

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

  adapt_wavelet ({f1, u}, (double[]){1e-3, 5e-2, 5e-2, 5e-2}, 
    minlevel = minlevel, maxlevel = maxlevel);

  solid_domain(cs, fs);
}

