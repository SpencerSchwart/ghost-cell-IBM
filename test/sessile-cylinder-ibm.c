/**
A 2D sessile droplet spreading on a cylindrical surface until it reaches an
equilibirum position, which is based on the contact angle. */


#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

/** 
Simulation variables */

const int maxlevel = 7;
const int minlevel = 4;

double df = 1;         // droplet diameter
double ds = 1;         // cylinder diameter
double theta0 = 0;     // contact angle
double Oh = 0.04;
double v0 = 0;         // initial volume

double t_end = 20;

FILE * outf; // out file

/** 
Boundary conditions */

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

u.n[bottom] = neumann(0);
pf[bottom] = dirichlet(0);
p[bottom] = dirichlet(0);

int main(int argc, char* argv[])
{
  size (4.);

  /**
  We shift the bottom boundary. */
  origin (0, -1.9);
  init_grid (1 << maxlevel);
  
  TOLERANCE = 1e-3;

  if (argc > 1)
    theta0 = atof(argv[1]);
  if (argc > 2)
    t_end = atof(argv[2]);
  if (argc > 3)
    Oh = atof(argv[3]);

  /**
  We use a constant fluid properties. */
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
    double angles[11] = {15,30,45,60,75,90.05,105,120,135,150,165};
    double tends[11] =  {30,15,15,15,15,15,15,15,20,20,20};
    for (int i = 0; i < 11; i++) {
      theta0 = angles[i];
      t_end = tends[i];

      if (theta0 >= 150)
        Oh = 0.16;

      mu1 = mu2 = Oh * sqrt(rho1 * f.sigma * df);

      f.wetting.theta_s = theta0;
      run();
    }
  }
}

void solid_domain (scalar cs, face vector fs)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (sq(x) + sq(y) - sq(ds/2.));
  boundary ({phi});
  fractions (phi, cs, fs);
  fractions_cleanup(cs,fs,1e-3);
}

event init (t = 0)
{
  /**
  We define the cylindrical wall and the initial (half)-circular
  interface. */
  solid_domain(cs, fs);
  fraction(f,  -(sq(x) + sq(y - (0.25*(df + ds))) - sq(df/2.)));
  fraction(ch, -(sq(x) + sq(y - (0.25*(df + ds))) - sq(df/2.)));

  v0 = real_volume(f, cs, fs);
}

event logfile (i++; t <= t_end)
{

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1. || f[] > cs[] - 1e-6)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-5)
    return true;

  /**
  Calculate volume, contact line position, and maximum capillary number.

  TODO: calculate angle and arc angle. */

  double vreal = 0, xcl = 100, ca = 0;
  foreach() {

    vreal += f[]*sq(Delta);

    if (on_interface(cs) && f[] > INT_TOL && f[] < cs[] - INT_TOL) {
      coord pint;
      if (is_contact_cell_bc(point, f, cs, fs, pint = &pint)) {
        pint = local_to_global_coord(point, pint);
        double angle = atan2(pint.y, pint.x);
        if (angle < xcl)
          xcl = angle;
      }
    }

    if (gc[] && cs[] == 1 && f[]) {
      double umag = sqrt(sq(u.x[]) + sq(u.y[]));
      if (umag > ca)
        ca = umag;
    }
  }
  double verror = i == 0? 0: (vreal - v0)/v0 * 100;
  ca *= mu1/f.sigma;

  if (i == 0) {
    char name[80];
    sprintf(name, "outs/out-%g", theta0);
    outf = fopen(name, "w");
  }

  fprintf(outf, "%d %g %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verror, xcl*180/pi, ca);
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
  output_facets_contact (f, ch, cs, fp);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (ch, kappa);
  foreach()
    if (cs[] < 1. || f[] > cs[] - 1e-6)
      kappa[] = nodata;

  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2*statsf(f).sum;

  sprintf(name, "results/results-%g", theta0);
  FILE * fp2 = fopen (name, "w");

  fprintf (fp2, "%d %g %.5g %.5g %.3g %g %g %g\n", N, theta0, R, R/sqrt(V/pi), s.stddev, 
                 v0, V, (V - v0)/v0 * 100);
  fclose(fp2);
}

/**
We use an intermediate field f1 to trick the adapt function into always
keeping the liquid interface at the maximum resolution. */

#if 1
event adapt (i++)
{
  scalar f1[], cs1[];
  foreach() {
    f1[] = f[];
    cs1[] = cs[];
  }
  adapt_wavelet ({cs1, f1, u}, (double[]){1e-3, 1e-3, 1e-3, 1e-3},
                           maxlevel = maxlevel, minlevel = minlevel);
  solid_domain(cs, fs);
}
#endif
