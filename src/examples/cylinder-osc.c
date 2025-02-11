#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
// #include "navier-stokes/double-projection.h" // why does this not work?
#include "view.h"

#define L0 29.257
#define D 1.
#define LEVEL 11

const int maxlevel = LEVEL;
const int minlevel = maxlevel - 5;
int Re;
const double U0 =  1.; // inlet velocity
const double t_end = 200;
const double tf_start = 100;
coord ci = {L0/4, L0/2}; // initial coordinates of cylinder
coord xc = {0, 0};
const double A = 0.2*D;
const double freq = 0.156;

face vector muv[];

u.n[left] = dirichlet ((U0));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top]   = neumann (0);
u.n[bottom] = neumann (0);

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (A*2*M_PI*freq*sin(2*M_PI*freq*(t)))

int main() {
  size(L0);
  init_grid (1 << (LEVEL - 2));
  mu = muv;
  TOLERANCE = 1.e-5; 
  CFL = 0.8;

  Re = 185;
  run();
}

event init (t = 0) {
  xc.y = -A*cos(2*M_PI*freq*t);
  solid (ibm, ibmf, sq(x - ci.x - xc.x) + sq(y - ci.y - xc.y) - sq(D/2));
  refine (ibm[] < 1 && ibm[] > 0 && level < maxlevel);

  foreach() {
    u.x[] = ibm[] * U0;
    u.y[] = 0;
  }
}

event moving_cylinder (i++) {
  xc.y = -A*cos(2*M_PI*freq*(t));
  solid (ibm, ibmf, sq(x - ci.x - xc.x) + sq(y - ci.y - xc.y) - sq(D/2));
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
  boundary ((scalar *) {muv});
}

scalar e[]; // pressure field diff. between timesteps
scalar p0[];
event logfile (i++, t <= t_end){

  double solidCells = 0, sgCells = 0;
  foreach(reduction(+:solidCells), reduction(+:sgCells)) {
    if (ibm[] == 0)
        solidCells += sq(Delta);
    if (ibm[] <= 0.5)
        sgCells += sq(Delta);
    e[] = p[] - p0[];
    p0[] = p[];
  }

  coord Fp;
  coord Fmu; 

  ibm_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));

  double vy = uibm_y (0,0,0);

  fprintf (stderr, "%d %g %d %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n",
           i, t, Re, mgpf.i, mgpf.nrelax, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, 
           CD, CL, Fp.x, Fmu.x, xc.y, vy, solidCells, sgCells);
}

event profile1 (t += 195.51) {
  char name[80];

  sprintf (name, "surface1-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach(serial) {
    if (ibm[] > 0 && ibm[] < 1) {
      coord interCell = {x,y}, midPoint, n;
      ibm_geometry(point, &midPoint, &n);

      foreach_dimension()
        midPoint.x = interCell.x + midPoint.x*Delta;

      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI) + 180;
      
      double cp = extrapolate_scalar (point, ibm, midPoint, n, p);
      double cp1 = image_pressure (point, p, midPoint);

      fprintf (fv5, "%g %g %g %g\n", xc.y, theta, cp, cp1);
    }
  }
  fflush (fv5);
  fclose (fv5);

  scalar omega[];
  vorticity (u, omega);

  sprintf (name, "%d-vort_final-%g.png", Re, t);
  view (fov = 2, tx = -0.26, ty = -0.50,
        width = 3000, height = 1500); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("ibm", "ibmf",lw = 5);
  save (name);
}

event profile2 (t += 197.12) {
  char name[80];

  sprintf (name, "surface2-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach(serial) {
    if (ibm[] > 0 && ibm[] < 1) {
      coord interCell = {x,y}, midPoint, n;
      double area = ibm_geometry(point, &midPoint, &n);

      foreach_dimension()
        midPoint.x = interCell.x + midPoint.x*Delta;

      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI) + 180;
      
      double cp = extrapolate_scalar (point, ibm, midPoint, n, p);
      double cp1 = image_pressure (point, p, midPoint);

      fprintf (fv5, "%g %g %g %g\n", xc.y, theta, cp, cp1);
    }
  }
  fflush (fv5);
  fclose (fv5);

  scalar omega[];
  vorticity (u, omega);

  sprintf (name, "%d-vort_final-%g.png", Re, t);
  view (fov = 2, tx = -0.26, ty = -0.50,
        width = 3000, height = 1500); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("ibm", "ibmf",lw = 5);
  save (name);
}

event adapt (i++) {
  scalar ibmsf[]; 
  foreach() {
    ibmsf[] = vertex_average (point, ibm);
  }
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,3e-4,3e-4},
		 maxlevel = LEVEL, minlevel = minlevel);
  event ("vof");
}


event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
