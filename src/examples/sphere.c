#include "grid/octree.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "view.h"
#include "lambda2.h"

#define L0 16.
#define D 1.

int maxlevel = 8;
int minlevel = 5;
int Re;
double U0 =  1.; // inlet velocity
double t_end = 150;
double tf_start = 20;
coord ci = {0,0,0};

face vector muv[];

u.n[left] = dirichlet ((U0));
u.t[left] = dirichlet (0);
u.r[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
u.r[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);
p[top]   = neumann (0);
pf[top]  = neumann (0);

u.n[bottom] = neumann (0);
p[bottom]   = neumann (0);
pf[bottom]  = neumann (0);

u.n[front] = neumann (0);
p[front]   = neumann (0);
pf[front]  = neumann (0);

u.n[back] = neumann (0);
p[back]   = neumann (0);
pf[back]  = neumann (0);

u_x_ibm_dirichlet (0);
u_y_ibm_dirichlet (0);
u_z_ibm_dirichlet (0);

int cdcount = 0;
double avgCD = 0, avgCL = 0;
int count = 0;

int main() {
  size(L0);
  init_grid (1 << (maxlevel - 2));
  origin ( -3.*D, -L0/2, -L0/2);
  mu = muv;

  TOLERANCE = 1.e-5; 
  CFL = 0.5;

  Re = 300;
  run();
}

event init (t = 0) {
  solid (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) + sq(z - ci.z) - sq(D/2));
  refine (ibm[] < 1 && level < maxlevel);
  solid (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) + sq(z - ci.z) - sq(D/2));

  foreach() {
    u.x[] = ibm[] * U0;
  }
}


#if 0
event moving_cylinder (i++) {
  solid (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) + sq(z - ci.z) - sq(D/2));
}
#endif

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

event logfile (i++; t <= t_end) {

  coord Fp, Fmu;
  ibm_force(p, u, muv, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x) / (0.5 * (U0) * sq(D/2.) * M_PI);
  double CL = (Fp.y + Fmu.y) / (0.5 * (U0) * sq(D/2.) * M_PI);
  double CS = (Fp.z + Fmu.z) / (0.5 * (U0) * sq(D/2.) * M_PI);
  
  fprintf (stderr, "%d %g %d %d %d %d %d %d %d %g %g %g\n",
          i, t, Re, mgpf.i, mgpf.nrelax, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax,
	  CD, CL, CS);
}

#if 1
event frequency (i++) {
  if (t >= tf_start && Re >= 290) {
    char name[80];
    sprintf (name, "freq-5-%d.dat", Re);
    FILE * fp = fopen (name, "a");
    foreach_point(5.0000,0.0000,0.0000) {
      fprintf (fp, "%d %g %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], u.z[], p[]);
    }
    fclose (fp);

    sprintf (name, "freq-7.5-%d.dat", Re);
    FILE * fp1 = fopen (name, "a");
    foreach_point(7.5000,0.0000,0.0000) {
      fprintf (fp1, "%d %g %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], u.z[], p[]);
    }
    fclose (fp1);
    count++;
  }
}
#endif

#if 1
scalar l2[];

event movies (t = tf_start; t += 0.25; t <= t_end) {
  if(Re == 300) {
  lambda2 (u, l2);
  scalar vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);

  char name[80];
  sprintf (name, "movie-%g.png", t);
  view (fov = 11.44, quat = {0.072072,0.245086,0.303106,0.918076},
	tx = -0.307321, ty = 0.22653, bg = {1,1,1},
	width = 802, height = 634);
  draw_vof ("ibm", "ibmf");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  save (name);
  }
}
#endif

event adapt (i++) {
  scalar ibmsf[];
  foreach()
    ibmsf[] = vertex_average (point, ibm);
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,2e-3,2e-2,2e-3},
		 maxlevel = maxlevel, minlevel = minlevel);
}

