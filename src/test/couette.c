/**
# Couette flow between rotating cylinders

We test embedded boundaries by solving the (Stokes) Couette flow
between two rotating cylinders. */

#include "grid/quadtree.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "view.h"

#define NAVIER 0
#define SLIP 0

double lambda = 0.2; // slip length

u.n[immersed] = dirichlet(0);

#if NAVIER
    u.t[immersed] = x*x + y*y > 0.14? navier_slip(lambda): dirichlet(distance(x, y));
#elif SLIP
    u.t[immersed] = x*x + y*y > 0.14? neumann(0): dirichlet(distance(x, y));
#else
    u.t[immersed] = x*x + y*y > 0.14? dirichlet(0): dirichlet(distance(x, y));
#endif

int main(int argc, char* argv[])
{
  /**
  Space and time are dimensionless. This is necessary to be able to
  use the 'mu = fm' trick. */
  mu = fm;

  size (1.25);
  DT = 1.;
  
  origin (-L0/2., -L0/2.);
  stokes = true;

  if (argc > 1)
    lambda = atof (argv[1]);

  for (N = 16; N <= 256; N *= 2)
    run();
}

scalar un[];

#define WIDTH 0.5

event init (t = 0) {
  /**
  The viscosity is unity. */
  
  /**
  The geometry of the embedded boundary is defined as two cylinders of
  radii 0.5 and 0.25. */

  solid (ibm, ibmf, difference (sq(0.5) - sq(x) - sq(y),
			     sq(0.25) - sq(x) - sq(y)));

  /**
  The outer cylinder is fixed and the inner cylinder is rotating with
  an angular velocity unity. */
  
  /**
  We initialize the reference velocity field. */
  
  foreach()
    un[] = u.y[];
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

/**
We compute error norms and display the angular velocity, pressure and
error fields using bview. */

#if SLIP
  #define powerlaw(r) (sq(0.25)/(sq(0.25) + sq(0.5))*(sq(0.5)/r + r))
#elif NAVIER
  #define powerlaw(r,lambda) (r*(0.5+lambda) + sq(0.5)*(lambda-0.5)/r)*(sq(0.25)/(sq(0.25)*(0.5+lambda)+sq(0.5)*(lambda-0.5)))
#else // NO-SLIP
  #define powerlaw(r) (r*(sq(0.5/r) - 1.)/(sq(0.5/0.25) - 1.))
#endif

event profile (t = end)
{
  scalar utheta[], e[], en[], es[];
  foreach() {
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    if (ibm[] > 0.) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
#if NAVIER
      e[]  = utheta[] - powerlaw(r,-lambda);
#else
      e[]  = utheta[] - powerlaw(r);
#endif
    }
    else
      e[] = p[] = utheta[] = nodata;
  }
  char name[80];

#if NAVIER
  sprintf(name, "navier-%g.dat", lambda);
#elif SLIP
  sprintf(name, "slip.dat");
#else
  sprintf(name, "noslip.dat");
#endif
  static FILE * fp = fopen(name, "w");

  norm n = normf(e);
  fprintf (fp, "%d %.7g %.7g %.7g %d %d %d %d %d %g\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, lambda);

  draw_vof ("ibm", "ibmf", filled = -1, fc = {1,1,1});
  squares ("utheta", spread = -1);
  save ("utheta.png");

  draw_vof ("ibm", "ibmf", filled = -1, fc = {1,1,1});
  squares ("p", spread = -1);
  save ("p.png");

  draw_vof ("ibm", "ibmf", filled = -1, fc = {1,1,1});
  squares ("e", spread = -1);
  save ("e.png");

  if (N == 32) {
#if NAVIER
    sprintf(name, "navier-prof-%g.dat", lambda);
#elif SLIP
    sprintf(name, "slip-prof.dat");
#else
    sprintf(name, "noslip-prof.dat");
#endif
    FILE * fpv = fopen(name, "w");
    foreach() {
      double theta = atan2(y, x), r = sqrt(x*x + y*y);
      if (fabs(utheta[]) < nodata)
          fprintf (fpv, "%g %g %g %g %g %g %g\n",
	           r, theta, u.x[], u.y[], p[], utheta[], e[]);
    }
    fclose(fpv);
  }

  fflush(fp);
  if (N == 256)
    fclose(fp);
}

