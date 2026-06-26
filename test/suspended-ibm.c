/**
File for simulating a 2D-axisymmetric suspended droplet on a flat, free-slip plate. 
The droplet is perturbed in a way in which its first mode is excited. */

#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-axi.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

/** 
Simulation variables */

const int maxlevel = 7;
const int minlevel = 3;

double df = 1;          // droplet diameter (5 mm)
double theta0 = 90.01;  // contact angle
double v0 = 0;          // initial volume

double omegac = 23.664; // capillary frequency (sigma/sqrt(rhol * D^3))
double gpert = 1.75;   // gravity perturbation to excite 1st mode (2 m/s^2)
double tpert = 0.56;    // time until perturbation ends (23.67 ms)

double tend = 40;
double tout = 0.25; 

FILE * outf; // out file

/** 
Boundary conditions */

// Free-slip wall
u.n[immersed] = dirichlet(0);
u.t[immersed] = neumann(0);

u.n[right] = neumann(0);
p  [right] = dirichlet(0);
pf [right] = dirichlet(0);

u.n[top] = neumann(0);

int main(int argc, char* argv[])
{
  size (2.);

  /**
  We shift the bottom boundary. */
  origin (-0.254, 0);
  init_grid (1 << minlevel);
 
  /**
  Set the multigrid solver tolerances. */
  TOLERANCE    = 1e-4;
  TOLERANCE_MU = 1e-3;

  /**
  We use fluid properties corresponding to water and air, normalized using the
  liquid density, surface tension, and the diameter of the droplet. */
  rho1 = 1;
  rho2 = 0.0012;
  mu1 = 1.69E-3; // Ohnesorge number
  mu2 = 3.04E-5;
  f.sigma = 1;

  f.wetting.theta_s = theta0;

  run();
}


void solid_domain (scalar cs, face vector fs)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = x;
  boundary ({phi});
  fractions (phi, cs, fs);
}


event init (t = 0)
{
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */

  scalar cs1[], f1[];
  face vector fs1[];

  int count = 0;
  astats st;
  do {
    solid_domain(cs1, fs1);
    fraction (f1, - (sq(x) + sq(y) - sq(df*0.5)));
    st = adapt_wavelet ({cs1, f1}, (double[]){1e-5, 1e-5}, maxlevel, minlevel);
    count++;
  } while ((st.nc || st.nf) && count < 40);

  solid_domain(cs, fs);
  event("update_metric");

  fraction (f,  - (sq(x) + sq(y) - sq(df/2.)));
  fraction (ch, - (sq(x) + sq(y) - sq(df/2.)));

  boundary(all);

  v0 = real_volume(f, cs, fs);
}


event acceleration (i++)
{
  if (t <= tpert) {
    face vector av = a;
    foreach_face(x)
      av.x[] -= (f[] + f[-1])*0.5*gpert; 
  }
}


event logfile (i++; t <= tend)
{

  /**
  Calculate volume, contact line position, and maximum capillary number. */

  double vreal = 0, xcl = 0, ca = 0;
  foreach(reduction(+:vreal) reduction(max:xcl) reduction(max:ca)) {

    vreal += f[]*sq(Delta);

    if (on_interface(cs) && f[] > INT_TOL && f[] < cs[] - INT_TOL) {
      coord pint;
      if (is_contact_cell_bc(point, f, cs, fs, pint = &pint)) {
        pint = local_to_global_coord(point, pint);
        if (pint.y > xcl)
          xcl = pint.y;
      }
    }

    if (gc[] && cs[] == 1) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]));
      if (umag > ca)
        ca = umag;
    }
  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  /**
  Measure the droplet height. */
  double hf = 0;
  foreach_boundary(bottom, reduction(max:hf)) {
    if (f[] > 1e-6 && f[] < 1-1e-6 && cs[] == 1) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double htemp = x + p.x*Delta;
      if (htemp > hf)
        hf = htemp;
    }
  }

  if (i == 0) {
    char name[80];
    sprintf(name, "outs/out-%g", theta0);
    outf = fopen(name, "w");
  }

  fprintf(outf, "%d %g %g %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verror, xcl, ca, hf);
  fflush(outf);
}

event end (t = end)
{
  fclose(outf);

  /**
  At the end, we output the equilibrium shape. */
  
  char name[80];
  sprintf (name, "shapes/shape-%g", theta0);
  FILE * fp = fopen (name, "w");
  output_facets (ch, fp);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;

  sprintf(name, "results/results-%g", theta0);
  FILE * fp2 = fopen (name, "w");

  fprintf (fp2, "%d %g %.5g %.3g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev, v0, V);
  fflush(fp2);
  fclose(fp2);
}


/**
We trick the adapt function to keep the liquid interface at the maximum level
of refinement. Note we do not do that with the solid interface. */
event adapt (i++) 
{
  scalar f1[];
  foreach()
    f1[] = ch[];

  adapt_wavelet ({cs, f1, u}, (double[]){1e-3, 1e-3, 5e-2, 5e-2}, 
    maxlevel = maxlevel, minlevel = minlevel);

  solid_domain (cs, fs);
}

