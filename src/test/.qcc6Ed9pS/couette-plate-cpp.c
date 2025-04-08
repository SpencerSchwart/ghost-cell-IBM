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
#line 1 "couette-plate.c"
//#include "navier-stokes/centered.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "view.h"


#define L0 10
int maxlevel = 6;

const double U0 = 5;
double Re = 50;
double t_end = 100;
double H = 10;
double offset = 2.45;
// coord vc = {0, 0};
const double delta = 10./pow(2,7);

// scalar vof[];
// face vector sf[];
face vector muv[];

u.n[left] = neumann (0);
//u.n[left] = dirichlet (1);

u.n[right] = dirichlet ( y < 1? 1: neumann (0));
//u.n[right] = neumann (0);
u.t[right] = neumann (0);
//p[right]   = dirichlet (0);
//pf[right]  = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (U0);
p[top] = neumann (0);
pf[top] = neumann (0);

// u.n[bottom] = dirichlet (0);
// u.t[bottom] = dirichlet (x <= 1? 0: neumann(0));
//u.t[bottom] = dirichlet (0);
 u.t[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

u_x_ibm_dirichlet(0)
u_y_ibm_dirichlet(0)

int main() {
  size(L0);
  init_grid (1 << (maxlevel));
  mu = muv;
  // TOLERANCE = 1.e-6;
  // TOLERANCE_MU = 1.e-5;
  // stokes = true;
  DT = 0.1;

  run();
}

double mm = 0.05;

event wall (i=0) {
  // solid (ibm, ibmf, y > (1.5) );
  //solid (ibm, ibmf, y > 2.55 && y < 7.45);
  solid (ibm, ibmf, y - offset);
}

/*
event ibm_advect (i++)
{
    solid (vof, sf, y > 1.5 + (t * 0.001));
}
*/



event properties (i++) {
  foreach_face() {
      //muv.x[] = fm.x[]*(U0)*(H)/(Re);  // viscosity must be zero in solid
        muv.x[] = ibmf.x[]*(U0)*(H)/(Re);  // viscosity must be zero in solid

  }
  boundary ((scalar *) {muv});
}

double v_profile (double x, double u, double offset) {
    return u * ((x-offset)/(H-offset));
}

event profile (t = t_end) {
  double delta = L0/pow(2,maxlevel);
  char name[80];

  double difmax = 1e-30;

  sprintf (name, "vprof");
  FILE * fv = fopen (name, "w");
  for (double i = 0; i <= L0; i += delta)
    foreach_point (L0/2, i) {
      double target = v_profile (y, U0, offset);
      double dif = target - u.x[];
      if (fabs(dif) > fabs(difmax) && (ibm[] != 0 || is_ghost_cell(point, ibm)))
        difmax = dif;
      if (ibm[] != 0 || is_ghost_cell(point, ibm))
        fprintf (fv, "%g %g %g %g %g\n", x, y, u.x[], target, dif);
    }
  fprintf (stderr, "Maximum difference = %g\n", difmax);
  fclose (fv);
}

event logfile (i++, t <= t_end) {
  fprintf (stderr, "%d %g %d %d %d %d %g\n", i, t,
           mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, TOLERANCE);
}

/*
event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 2);
}
*/



event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%g\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fprintf (fp, "%d\t%d\t%d\t%d\n", mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  fflush (fp);
  return 1;
}


#endif
