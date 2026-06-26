/**
A 2D sessile droplet oscillating on a flat, horizontal plate until it reaches an
equilibirum position, which is based on the contact angle. The oscillation occurs
from low Ohnesorge number and initalizing the droplet at 90 degrees instead of
its equilibirum contact angle (60 or 120). */

#define GCV 1e-4

#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

/** 
Simulation variables */

int maxlevel = 6;
int minlevel = 4;

double df = 1;         // droplet diameter
double theta0 = 0;     // contact angle
double Oh = 0.005;     // Ohnesorge number
double v0 = 0;         // initial volume

double doff = 0.2;     // solid volume fraction in interfacial cells

double t_end = 40;

FILE * outf; // out file

/** 
Boundary conditions */

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

u.n[top] = neumann(0);
p[top] = dirichlet(0);
pf[top] = dirichlet(0);

int main(int argc, char* argv[])
{
  size (2.);

  if (argc > 1)
    theta0 = atof(argv[1]);
  if (argc > 2)
    doff = atof(argv[2]);
  if (argc > 3)
    t_end = atof(argv[3]);
  if (argc > 4)
    maxlevel = atoi(argv[4]);

  /**
  We shift the bottom boundary. */
  double hmin = 2./(1 << maxlevel);
  origin (0, -0.25-(doff*hmin));
  init_grid (1 << maxlevel);
  
  TOLERANCE = 1e-3;

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
     double angles[2] = {110,70};
     double doffs[4] = {0.2, 0.4, 0.6, 0.8};
     for (int i = 0; i < 2; i++) {
       theta0 = angles[i];
       f.wetting.theta_s = theta0;
       for (int j = 0; j < 4; j++) {
         doff = doffs[j];
         origin (0, -0.25-(doff*hmin));
         run();
       }
     }
  }
}

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
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  solid_domain(cs, fs);
  fraction (f, - (sq(x) + sq(y) - sq(df/2.)));
  fraction (ch, - (sq(x) + sq(y) - sq(df/2.)));

  v0 = real_volume(f, cs, fs);
}

event logfile (i++; t <= t_end)
{

  /**
  This block of code outputs a snapshot of the domain used for the paper. */
#if 0
  if (i == 5) {
    view(tx = -0.5, ty = -0.35, fov = 20, width = 2000, height = 2000);
    cells();
    draw_vof("cs", "fs", filled = -1, lw = 2, lc={0,0,0}, fc={0.796875, 0.796875, 0.796875});
    draw_vof("ch", filled = 1, lw = 4, lc={0,0,0}, fc={1, 0.62891, 0.62891});
    save("domain.png");
  }
#endif

  /**
  Calculate volume, contact line position, and maximum capillary number. */

  double vreal = 0, xcl = 0, ca = 0;
  foreach() {

    vreal += f[]*sq(Delta);

    if (on_interface(cs) && f[] > INT_TOL && f[] < cs[] - INT_TOL) {
      coord pint;
      if (is_contact_cell_bc(point, f, cs, fs, pint = &pint)) {
        pint = local_to_global_coord(point, pint);
        if (pint.x > xcl)
          xcl = pint.x;
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
  Calculate the height of the droplet. */
  double hf = 0;
  foreach_boundary(left) {
    if (f[] > 0 && f[] < 1 && cs[] == 1) {
      coord n = interface_normal(point, f), mp;
      double alpha = plane_alpha(f[], n);
      plane_area_center (n, alpha, &mp);
      mp = local_to_global_coord(point, mp);
      if (mp.y > hf)
        hf = mp.y;
    }
  }

  if (i == 0) {
    char name[80];
    sprintf(name, "outs/%g-out-%g-%g-%d", GCV, doff, theta0, maxlevel);
    outf = fopen(name, "w");
  }

  fprintf(outf, "%d %g %g %g %g %g %g %g %g %g\n", i, t, theta0, doff, v0, vreal, verror, xcl, ca, hf);
  fflush(outf);
}

/**
Given the radius of curvature $R$ and the volume of the droplet $V$,
this function returns the equivalent contact angle $\theta$ verifying
the equilibrium solution
$$
V = R^2(\theta - \sin\theta\cos\theta)
$$
*/

double equivalent_contact_angle (double R, double V)
{
  double x0 = 0., x1 = pi;
  while (x1 - x0 > 1e-4) {
    double x = (x1 + x0)/2.;
    double f = V - sq(R)*(x - sin(x)*cos(x));
    if (f > 0.)
      x0 = x;
    else
      x1 = x;
  }
  return (x0 + x1)/2.;
}

event end (t = end)
{
  fclose(outf);

  /**
  At the end, we output the equilibrium shape. */
  
  char name[80];
  sprintf (name, "shapes/%g-shape-%g-%g-%d", GCV, doff, theta0, maxlevel);
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
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;

  sprintf(name, "results/%g-results-%g-%g-%d", GCV, doff, theta0, maxlevel);
  FILE * fp2 = fopen (name, "w");

  fprintf (fp2, "%d %g %.5g %.3g %.4g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi, v0, V);
  fclose(fp2);
}

/**
We use an intermediate field f1 to trick the adapt function into always
keeping the liquid interface at the maximum resolution. */

event adapt (i++)
{
  scalar f1[];
  foreach()
    f1[] = ch[];
  adapt_wavelet ({f1, u,}, (double[]){1e-3, 1e-2, 1e-2},
                           maxlevel = maxlevel, minlevel = minlevel);
  solid_domain(cs, fs);
}
