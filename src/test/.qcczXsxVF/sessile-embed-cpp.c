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
#line 1 "sessile-embed.c"
/**
# Sessile drop on an embedded boundary

This test case is very close to the [sessile
drop](/src/test/sessile.c) test but uses an embedded boundary rather
than the boundary of the domain.

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1.6:1.6]
set yrange [0:]
plot 'out' w l, '' u (-$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/
#define LIMIT 1e10
#define on_interface(a) (a[] > 1e-6 && a[] < 1-1e-6)
#define distance(a,b) sqrt(sq(a) + sq(b))

#include "embed.h"
#include "../my-centered.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-embed.h"

double theta0;

bool is_triple_point (Point point, coord nf, coord ns)
{
    if (!(on_interface(cs)) || !(on_interface(f)))
        return false;
    if (ns.x == 0 || ns.y == 0 || nf.x == 0 || nf.y == 0)
        return false;

    double alphas = plane_alpha (cs[], ns);

    double alphaf = plane_alpha (f[], nf);

    double intercept = ((alphas/ns.y) - (alphaf/nf.y)) /
                       ((ns.x/ns.y) - (nf.x/nf.y));
    
    return fabs(intercept) <= 0.5;

}

double get_contact_angle (scalar f, scalar ibm)
{
    vector nf[], ns[];
    scalar alphaf[], alphas[];

    reconstruction (f, nf, alphaf);
    reconstruction (ibm, ns, alphas);

    double theta = 0;
    int count = 0;
    foreach() {
        if (on_interface(ibm) && on_interface(f)) {
            coord nf_temp = {nf.x[], nf.y[]}, ns_temp = {ns.x[], ns.y[]};
            if (is_triple_point (point, nf_temp, ns_temp)) {
                double num = nf.x[]*ns.x[] + nf.y[]*ns.y[];
                double den = distance(nf.x[], nf.y[]) * distance(ns.x[], ns.y[]);
                theta += acos (num/den);
                count++;
            }
        }
    }

    if (count > 0)
        return (theta / count)*180./pi;
    else
        return 0;
}

int main()
{
  size (2.);

  /**
  We shift the bottom boundary. */

  origin (0, - 0.26);
  
  /**
  We use a constant viscosity. */

  mu1 = mu2 = 0.1;
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  We vary the contact_angle. */

  for (theta0 = 15; theta0 <= 165; theta0 += 15) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
}

double v0 = 0;

event init (t = 0)
{

  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = y;
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));

  foreach()
    v0 += f[]*sq(Delta);

}

event logfile (i++; t <= 20)
{

  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;

  double vreal = 0, vreal2 = 0;
  foreach() {
    vreal += cs[]? f[]*sq(Delta): 0;
    vreal2 += f[]*dv();
  }
  
  double thetar = get_contact_angle(f, cs);

  fprintf(stderr, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, vreal2, thetar);
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
    if (cs[] < 1.)
      kappa[] = nodata;
  
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;

  double thetar = get_contact_angle(f, cs);

  fprintf (stdout, "%d %g %.5g %.3g %.4g %g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi), thetar;
}

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

The accuracy is not as good (yet) as that of the [sessile
test](/src/test/sessile.c).

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
contact angle" with the imposed contact angle.

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

* [Sessile drop on an inclined boundary](sessile-inclined.c)
* [Sessile drop on a domain boundary](/src/test/sessile.c)
*/


#endif
