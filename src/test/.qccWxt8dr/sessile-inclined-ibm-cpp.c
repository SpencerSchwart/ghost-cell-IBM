@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "sessile-inclined-ibm.c"

/**
# Sessile drop on an embedded boundary inclined at 45 degrees

This is a less trivial version of [this test case](sessile.c).

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 480,480
set size ratio -1
unset key
unset xtics
unset ytics
unset border
plot 'out' w l, x + 0.97 with filledcurves x1 lc 'grey'
set term pop
~~~
*/
#undef dv()
#define dv() (sq(Delta) * ibm[])
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm.h"

face vector acc[];

double theta0;

u_x_ibm_dirichlet(0)
u_y_ibm_dirichlet(0)

double real_volume (scalar c)
{
    vector nf0[], ns0[];
    scalar alphaf0[], alphas0[];
    reconstruction (c, nf0, alphaf0);
    reconstruction (ibm, ns0, alphas0);

    double sum = 0;
    foreach() {
        if (on_interface(ibm) && on_interface(c)) {
            coord nft = {nf0.x[], nf0.y[]}, nst = {ns0.x[], ns0.y[]};
            coord lhs = {-0.5, -0.5}, rhs = {0.5, 0.5};
            sum += immersed_area(c[], nft, alphaf0[], nst, alphas0[], lhs, rhs, 0)*(sq(Delta));
        }
        else
            sum += ibm[]*c[]*sq(Delta);
    }

    return sum;
}

int main()
{
  size (2.);
  origin (-1, 0);
  
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */
 
#if 0
  for (theta0 = 15; theta0 <= 165; theta0 += 15) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
#endif
  theta0 = 15;
  const scalar c[] = theta0*pi/180;
  contact_angle = c;
  run();
}

double init_volume = 0, v0 = 0;
event init (t = 0)
{
  /**
  We define the inclined wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - x - 0.97);
  boundary ({phi});
  fractions (phi, ibm, ibmf);
  fraction (f, - (sq(x - 0) + sq(y - 1.) - sq(0.25)));

  init_volume = 0;
  foreach() 
    init_volume += f[]*dv();
  v0 = real_volume(f);
}

event logfile (i++; t <= 20)
{
  foreach_face()
    acc.x[] = a.x[];
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;
  double vreal = real_volume(f), vreal2 = 0;
  foreach()
    vreal2 += cr[]*sq(Delta);

  fprintf(stderr, "%d %g %g %g %g\n", i, t, v0, vreal, vreal2);
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
  double final_volume = 0;
  foreach()
    final_volume += f[]*dv();
    //final_volume += f[]*(ibm[] != 0)*sq(Delta);
  /**
  At the end, we output the equilibrium shape. */
  
  char name[80];
  sprintf (name, "shape-%g", theta0);
  FILE * fp = fopen (name, "w");
  output_facets (f, fp);

  /**
  We compute the curvature only in full cells. */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (ibm[] < 1.)
      kappa[] = nodata;
    
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.4g %g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi, init_volume, final_volume);
}

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

The results are not very good, especially for small contact angles.

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  'log' u 2:3 pt 7 t 'numerical'
~~~

Another way to display the same result is to compare the "apparent
contact angle" with the imposed contact angle. It is clear that
contact angles smaller than $\approx 45$ degrees or larger than
$\approx 135$ degrees cannot be imposed.

~~~gnuplot
reset
set xlabel 'Imposed contact angle (degrees)'
set ylabel 'Apparent contact angle (degrees)'
unset key
set xtics 15,15,165
set ytics 15,15,165
plot 'log' u 2:5 pt 7, x
~~~

## See also

* [Sessile drop on a horizontal embedded boundary](sessile.c)
*/

#endif
