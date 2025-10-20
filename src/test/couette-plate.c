//#include "navier-stokes/centered.h"
#include "../ibm-gcm.h"
#define CCM 0
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "view.h"

#define L0 10
const int l0 = 10;
const int maxlevel = 6;
//const double delta = l0/pow(2,maxlevel);

const double U0 = 5;
double Re = 50;
double t_end = 100;
double H = 10;
double offset = 2.4;

face vector muv[];

u.n[left] = neumann (0);
u.t[left] = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (U0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.t[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

//u.t[immersed] = neumann(0);
//u.t[immersed] = navier_slip(0.2);
u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

int main() {
  size(L0);
  init_grid (1 << (maxlevel));
  mu = muv;
  // TOLERANCE = 1.e-6;
  // TOLERANCE_MU = 1.e-5;
  stokes = true;
  DT = 0.1;

  run();
}

event wall (i = 0) {
  solid (ibm, ibmf, y - offset);
}

event properties (i++) {
    foreach_face() {
        muv.x[] = fm.x[]*(U0)*(H)/(Re);  // viscosity must be zero in solid
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

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%g\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fprintf (fp, "%d\t%d\t%d\t%d\n", mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  fflush (fp);
  return 1;
}

