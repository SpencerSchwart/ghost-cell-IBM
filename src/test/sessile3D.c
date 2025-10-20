/**
# 3D Sessile drop

This is the 3D equivalent of the [2D](sessile.c) test case.

The volume of a [spherical
cap](https://en.wikipedia.org/wiki/Spherical_cap) of radius $R$ and
(contact) angle $\theta$ is
$$
V = \frac{\pi}{3}R^3(2+\cos\theta)(1-\cos\theta)^2
$$
or equivalently
$$
\frac{R}{R_0} = \left(\frac{1}{4}(2+\cos\theta)(1-\cos\theta)^2\right)^{-1/3}
$$
with $R_0$ the equivalent radius of the droplet
$$
R_0 = \left(\frac{3V}{4\pi}\right)^{1/3}
$$
To test this relation, a drop is initialised as a half-sphere
(i.e. the initial contact angle is 90$^\circ$) and the contact angle
is varied between 30$^\circ$ and 150$^\circ$. The drop oscillates and
eventually relaxes to its equilibrium position. The curvature along
the interface is close to constant.

![Relaxation toward a $30^\circ$ contact angle.](sessile3D/movie.mp4)

Note that shallower angles are [not accessible
yet](/src/contact.h). */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "vof.h"
#include "navier-stokes/perfs.h"

#include "tension.h"
#include "view.h"

#define L0 2.
#define D0 1.


double theta0 = 90;
vector h[];
h.t[bottom] = contact_angle (theta0*pi/180.);
h.r[bottom] = contact_angle (theta0*pi/180.);

u.n[bottom] = dirichlet(0);
u.t[bottom] = neumann(0);
u.r[bottom] = neumann(0);

u.n[left] = dirichlet(0);
u.t[left] = neumann(0);
u.r[left] = neumann(0);

u.n[back] = dirichlet(0);
u.t[back] = neumann(0);
u.r[left] = neumann(0);

f[left] = neumann(0);
f[back] = neumann(0);

p[left] = neumann(0);
p[back] = neumann(0);

pf[left] = neumann(0);
pf[back] = neumann(0);

p[top]  = dirichlet(0);
pf[top] = dirichlet(0);

int max_level = 8;
int min_level = 3;
double L = 2.;
double t_end = 100;
double amp = 0.e-2;
double dia = 5.e-3;
double gpert = 4.9;

double femax = 0.01;
double uemax = 0.001;

double rhol  = 1000;
double rhog  = 1.2;
double mul   = 1.e-3;
double mug   = 1.8e-5;
double sigma = 0.07;

double omegac = 0.5367;
double time_pert = 1.0434;

int main()
{
  size (L);
  init_grid (1 << (max_level - 3));
  
  TOLERANCE = 1e-5;
  rho1 = rhol/rhol;
  rho2 = rhog/rhol;
  mu1 = mul/sqrt(rhol*dia*sigma);
  mu2 = mug/sqrt(rhol*dia*sigma);

  f.sigma = sigma/sigma;

  f.height = h;

  run();
}

/**
The initial drop is a quarter of a sphere. */

event init (t = 0)
{
  scalar f1[];
  fraction (f1, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
  refine (f1[] > 0 && f1[] < 1 && level < max_level);
  fraction (f, - (sq(x) + sq(y) + sq(z) - sq(D0/2.)));
}

/**
We log statistics on the maximum velocity, curvature and volume. If
the standard deviation of curvature falls below $10^{-2}$, we assume
that the steady shape is reached and we stop the calculation. */

event logfile (i ++; t <= t_end)
{
#if 0
  scalar kappa[];
  cstats cs = curvature (f, kappa);
  foreach()
    if (f[] <= 1e-3 || f[] >= 1. - 1e-3)
      kappa[] = nodata;
  stats s = statsf (kappa);
#endif

  scalar pos[];
  position (f, pos, {0,1,0});
  double hmax = statsf(pos).max;
  position (f, pos, {1,0,0});
  double rmax = statsf(pos).max;

#if 0
  fprintf (fout, "%d %g %g %g %g %g %g %g %g %d %d %d %d %g %g\n", i, t, theta0, normf(u.x).max,
	   s.min, s.sum/s.volume, s.max, s.stddev, statsf(f).sum,
	   cs.h, cs.f, cs.a, cs.c, hmax, rmax);
#else
  fprintf (fout, "%d %g %g %g %g\n", i, t, theta0, hmax, rmax);
#endif
  fflush (fout);
}

event acceleration (i++)
{
  foreach_face(y)
    a.y[] += (t > time_pert? 0:
               -gpert*rhol*sq(dia)/sigma);
}

#if 0
event snapshots (i += 10)
{
  scalar kappa[];
  curvature (f, kappa);
  p.nodump = false;
  dump (buffered = true);
}
#endif

/**
We make a movie of the relaxing interface for $\theta = 30^\circ$. We
use symmetries since only a quarter of the drop is simulated. */

event movie (t += 0.1; t <= t_end)
{
char name[80];
    view (fov = 26.6776, quat = {0.474458,0.144142,0.234923,0.836017},
	  tx = -0.0137556, ty = -0.00718937, bg = {1,1,1},
	  width = 758, height = 552);
    draw_vof ("f");
    draw_vof ("f", edges = true);
    cells (lc = {1,0,0});
    mirror (n = {1,0,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
    }
    mirror (n = {0,1,0}) {
      draw_vof ("f");
      draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
      mirror (n = {1,0,0}) {
	draw_vof ("f");
	draw_vof ("f", edges = true);
	cells (lc = {1,0,0});
      }
    }
    save ("movie.mp4");

  view (quat = {-0.280, 0.155, 0.024, 0.947}, fov = 30, near = 0.01, far = 1000, bg = {1,1,1},
      tx = 0, ty = 0.299, tz = -2.748, width = 700, height = 550);
  cells (n = {0,1,0});
  draw_vof (c = "f", fc = {0.647,0.114,0.176}, lw = 2);
  draw_vof (c= "f", edges = true);
  sprintf(name, "%g-f-%g.png", theta0, t);
  save (name);
}

/**
We use refinement based on a smooth version of the volume
fraction. This guarantees constant refinement around the interface,
which seems to be necessary to reach balance (this should be
improved). Note however that this is specific to this test case and
should not generally be used in "production" runs, which should work
fine with the default criterion. */

event adapt (i++) {
  scalar f1[];
  foreach()
    f1[] = f[];
  adapt_wavelet ({f1, u}, (double[]){1e-3, 3e-2, 3e-2, 3e-2}, 
    minlevel = min_level, maxlevel = max_level);
}
