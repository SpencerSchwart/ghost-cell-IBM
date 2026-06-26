#include "grid/octree.h"
#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "view.h"
#include "lambda2.h"

#include "navier-stokes/perfs.h"

double l0 = 16;
double ds = 1;

int maxlevel = 9;
int minlevel = 5;

int Re = 0;
double u0 =  1.; // inlet velocity

double tend = 100.01;
double tfstart = 50;
double tout = 0.25;

face vector muv[]; // viscosity

FILE * fout;

/**
Boundary conditions. */

/* Inlet */
u.n[left] = dirichlet (u0);
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

/* Outlet */
u.n[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/* No-Slip Wall */
u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

int main(int argc, char* argv[]) {

  if (argc > 1)
    Re = atoi(argv[1]);
  if (argc > 2)
    tend = atof(argv[2]);
  if (argc > 3)
    maxlevel = atoi(argv[3]);

  size(l0);
  init_grid (1 << (minlevel));
  origin ( -l0/4, -l0/2, -l0/2);

  mu = muv;

  /**
  Low tolerance to help get less noisy $C_p$ distribution and pressure drag. The
  viscosity solver's tolerance is kept to the default value. */

  TOLERANCE    = 1e-5; 
  TOLERANCE_MU = 1e-3;

  CFL = 0.8;

  if (Re) {
    run();
  }
  else {
    tend = 50.001;
    Re = 100;
    run();

    tend = 100.001;
    Re = 250;
    run();

    tend = 150.001;
    Re = 300;
    run();
  }
}

/**
Since the geometry is trivial, we can use the exact normal and boundary intercept
values for imposing the no-slip BC on the solid's surface to improve the accuracy.

TODO: how much does this effect the solution? */

coord solid_normal (coord p) 
{
  coord n = {p.x, p.y, p.z};
  normalize_sum(&n);
  return n;
}

coord solid_bi (coord p)
{
  double d = magnitude_coord(p);
  coord bi = {0};
  foreach_dimension()
      bi.x = (0.5*ds*p.x)/d;
  return bi;
}

scalar cf_avg[];
scalar cp_avg[];

int cf_counter = 0;

event init (t = 0) 
{
  
  /**
  Iterative refinement to initalize the rigid sphere. Calling the update_metric event
  after seems to avoid crashing in the stability event caused by uf != 0 when fm == 0. */

  scalar cs1[];
  face vector fs1[];

  astats ad;
  int count = 0;
  do {
      solid (cs1, fs1, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
      ad = adapt_wavelet ({cs1}, (double[]){1e-5}, minlevel = minlevel, maxlevel = maxlevel);
      count++;
  } while ((ad.nc || ad.nf) && count < 20);
  
  solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
  event("update_metric");

  boundary ({cs});

  /**
  Initialize the velocity field and the average C_f and C_p distribution fields. */

  foreach() {
    cf_avg[] = cp_avg[] = 0;
    u.x[] = cs[]*u0;
  }

  foreach_face()
    muv.x[] = fm.x[]*(u0)*(ds)/(Re);
   boundary ((scalar *) {muv});

  restriction(all);
  boundary(all);

  cf_counter = 0;

  /**
  Set the exact normal and boundary intercept functions. */

  cs.ibm.normal = solid_normal;
  //cs.ibm.bi = solid_bi;

  char name[80];
  sprintf(name, "outs/%d-out-%d", maxlevel, Re);
  fout = fopen(name, "w");
}

/**
We overload the properties event to set the viscosity based on the Reynolds number. */

event properties (i++) 
{
  foreach_face()
    muv.x[] = fm.x[]*(u0)*(ds)/(Re);
   boundary ((scalar *) {muv});
}

/**
We calculate the drag, lift, and side force coefficients. */

event logfile (i++; t <= tend)
{

  vector ue[];
  foreach() {
    foreach_dimension()
      ue.x[] = nodata;
    if (cs[] > 0 && cs[] < 1) {
      coord n = facet_normal (point, cs, fs), midp;
      double alpha = plane_alpha (cs[], n);
      plane_area_center (n, alpha, &midp);
      midp = local_to_global_coord(point, midp);

      ue.x[] = interpolate_linear(point, u.x, midp.x, midp.y, midp.z);
      ue.y[] = interpolate_linear(point, u.y, midp.x, midp.y, midp.z);
      ue.z[] = interpolate_linear(point, u.z, midp.x, midp.y, midp.z);
    }
  }

  norm2 nue = normf2 (ue.x);
  norm2 nve = normf2 (ue.y);
  norm2 nwe = normf2 (ue.z);

  coord Fp, Fmu;
  ibm_force(p, u, muv, &Fp, &Fmu);

  // pressure force
  double CDp = (Fp.x) / (0.5 * (u0) * (0.25*sq(ds)*pi));
  double CLp = (Fp.y) / (0.5 * (u0) * (0.25*sq(ds)*pi));
  double CSp = (Fp.z) / (0.5 * (u0) * (0.25*sq(ds)*pi));

  // viscous/friction force
  double CDf = (Fmu.x) / (0.5 * (u0) * (0.25*sq(ds)*pi));
  double CLf = (Fmu.y) / (0.5 * (u0) * (0.25*sq(ds)*pi));
  double CSf = (Fmu.z) / (0.5 * (u0) * (0.25*sq(ds)*pi));

  norm2 nnorm = normf2(neg);
  stats nstat = statsf(neg);

  fprintf (fout, "%d %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
    i, t, Re, CDp, CLp, CSp, CDf, CLf, CSf, nue.avg, nue.rms, nue.max,
    nve.avg, nve.rms, nve.max, nwe.avg, nwe.rms, nwe.max, 
    nnorm.avg, nnorm.max, nstat.stddev);
  fflush(fout);
}

/**
For Re >= 300, vortex shedding occurs. We want to capture this frequency
to calculate the Strohaul number. After an initial time tfstart, velocity
is measured at a probe 5*D down stream. Preforming an FFT on this data
yields the frequency. */

int count = 0;
event frequency (i++) 
{
  static const coord probe = {5.00000, 0.00000, 0.00000};

  if (t >= tfstart && Re >= 300) {
    double xp = -nodata, yp = -nodata, zp = -nodata, up = -nodata, vp = -nodata, wp = -nodata;

    foreach_point(probe.x * ds, probe.y, probe.z, serial) {
        xp = x;
        yp = y;
        zp = z;
        up = interpolate_linear (point, u.x, probe.x * ds, probe.y, probe.z);
        vp = interpolate_linear (point, u.y, probe.x * ds, probe.y, probe.z);
        wp = interpolate_linear (point, u.z, probe.x * ds, probe.y, probe.z);
    }
    
    if (xp > -nodata) {
        char name[80];
        sprintf (name, "data/%d-freq-%d.dat", maxlevel, Re);
        FILE * fp = fopen (name, "a");
        fprintf (fp, "%d %g %g %g %g %g %g %g\n", count, t, xp, yp, zp, up, vp, wp);
        fclose (fp);
        count++;
    }
  }
}

/**
We are interested in outputting the distribution of $C_p$ and $C_f$ in the x-y plane.
For unsteady Re, we should average it in time, though we do it for all Re anyways. */

event surface_data_avg (t = tfstart; i += 5)
{
  scalar cf[];
  skin_friction (u, mu, cf);
  
  foreach() {
    if (cs[] > 0 && cs[] < 1) {
        coord midPoint, n, b, cc = {x,y,z};
        ibm_geometry (point, &b, &n);
  
        foreach_dimension()
          midPoint.x = cc.x + b.x*Delta;

        cp_avg[] += extrapolate_scalar (point, cs, midPoint, n, p) / (0.5 * sq(u0));
        cf_avg[] += cf[] / (0.5 * sq(u0));
    }
  }
  cf_counter++;
}

/**
Here we output the final and average distribution along the surface in the x-y plane. */

event surface_data (t = end)
{
  char name[80];
  coord pc;
  coord xy_box[2] = { {-1, -1, 0}, {1, 1, 0} }; // x-y plane
  coord n = {512,512};

  sprintf (name, "data/%d-%d-xysurface-%d", pid(), maxlevel, Re);
  FILE * fpxy = fopen (name, "w");

  sprintf (name, "data/%d-%d-xysurface-avg-%d", pid(), maxlevel, Re);
  FILE * fpxy_avg = fopen (name, "w");

  scalar cf[];
  skin_friction (u, mu, cf);
  
  foreach_region (pc, xy_box, n, serial) {

    if (cs[] > 0 && cs[] < 1) {

      coord midPoint, n, b, cc = {x,y,z};
      double area = ibm_geometry (point, &b, &n);

      foreach_dimension()
        midPoint.x = cc.x + b.x*Delta;

      double theta = atan2 (y,x) * (180/M_PI) + 180;
      double cp = extrapolate_scalar (point, cs, midPoint, n, p) / (0.5 * sq(u0));
      double cfi = cf[] / (0.5 * sq(u0));

      fprintf (fpxy, "%g %g %g %g %g %g %g\n", 
                      pc.x, pc.y, pc.z, theta, area, cp, cfi);

      fprintf (fpxy_avg, "%g %g %g %g %g %g %g %g %g %g\n", 
                      midPoint.x, midPoint.y, midPoint.z, pc.x, pc.y, pc.z, theta, area, 
                      cp_avg[]/((double)cf_counter + SEPS), cf_avg[]/((double)cf_counter + SEPS));
    }
  }
  
  fclose(fpxy);
  fclose(fpxy_avg);

  sprintf(name, "data/%d-fields-%d", maxlevel, Re);
  FILE * fpf = fopen(name, "w");
  output_field((scalar *){u.x, u.y, p, cs}, fp = fpf, n = 1 << maxlevel);
  fclose(fpf);

  sprintf(name, "data/%d-interface", maxlevel);
  FILE * fpi = fopen(name, "w");
  output_facets(cs, fp = fpi, s = fs);
  fclose(fpi);

  scalar umag[];
  foreach()
    umag[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
  sprintf(name, "imgs/%d-umag-%d", maxlevel, Re);

  view (fov = 2.5, width = 2000, height = 500);
  squares("umag", min=-1.5, max=1.5);
  draw_vof("cs", "fs", fc={0.8, 0.8, 0.8});
  draw_vof ("cs", edges = true, lw = 0.5);
  cells(lw = 0.5);
  save(name);

  fclose(fout);
}

/**
Generate a movie with the lambda2 isosurfaces whose color represents the y-z vorticity. */

scalar l2[];

event movies (t = tfstart; t += tout) 
{
  if (Re == 300) {
    lambda2 (u, l2);
    scalar vyz[];
    foreach()
      if (gc[])
        vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  
    char name[80];
    sprintf (name, "imgs/%d-lambda2-%.5f.png", Re, t);
    
    view (fov = 11.44, quat = {0.072072,0.245086,0.303106,0.918076},
  	      tx = -0.307321, ty = 0.22653, bg = {1,1,1}, width = 1000, height = 800);
    draw_vof ("cs", "fs");
    isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1, linear = true, map = cool_warm);
    save (name);
  }
}

event dumps (t += tout)
{
  p.nodump = false;
  scalar * dumplist = {u,p,cs};
  char name[80];
  sprintf(name, "dumps/%d-dump-%g", Re, t);
  dump (file = name, list = dumplist);
}

/**
We adapt based on both u and cs so that the solid surface is at the maximum level of
refinement. This facilitates a less noisy surface force calculation. To avoid AMR crashes, 
we reinitialize the solid surface after every adaption. */

event adapt (i++) 
{
  scalar cs1[];
  foreach()
    cs1[] = cs[];
  adapt_wavelet ({cs1, u}, (double[]){1e-3, 5e-3, 5e-3, 5e-3},
		 maxlevel = maxlevel, minlevel = minlevel);

  solid (cs, fs, sq(x) + sq(y) + sq(z) - sq(0.5*ds));
}

