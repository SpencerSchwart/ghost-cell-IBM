#define CA 1
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm3D.h"
#include "view.h"

#define LEVEL 8
#define RHOR 80

int j = 3, k = 0;
double uemax = 1;
double t_end = 2.5 [0,1];

const double H = 1 [1];  // gas nozzle 
const double R = 3*H;  // liquid inlet
const double RHO1 = 1;
const double SIGMA = 1.;

double MUl = 5.e-3;
double MUg = 2.5e-2;
double Ul = 0.5;
double Ug = 160;
double theta0 = 70;

u.n[left] = dirichlet (fabs(y) < R? Ul:0);
u.t[left] = dirichlet (0);
p[left] = neumann (0);
pf[left] = neumann (0);

u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

u.n[top] = dirichlet (x > 0 && x <= H? -Ug:0);
u.t[top] = dirichlet (0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = dirichlet (x > 0 && x <= H? Ug:0);
u.t[top] = dirichlet (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (0)

f[left] = dirichlet (fabs(y) < R? 1: 0);

// Define solid/nozzle shape
static void solid_domain (scalar c, face vector cf) {
  vertex scalar phi[];
  foreach_vertex() {
    if ((x >= H && y >= R) || (x <= 0 && y >= R) || 
	(x <= 0 && y <= -R) || (x >= H && y <= -R))
      phi[] = -1.;
    else 
      phi[] = 1.;
  }
  boundary ({phi});
  fractions (phi, c, cf);
  fractions_cleanup (c, cf);
}

// Define which domain parts are liquid vs gas
static void lq_domain (scalar c) {
  vertex scalar phi[];
  foreach_vertex () {
   if (x < -H)
      phi[] = 1.;
    else
      phi[] = -1.;
  }
  boundary ({phi});
  fractions (phi, c);
}


int main (int argc, char * argv[]) {
  size (25*H);
  init_grid (1 << (LEVEL - 2));
  origin (-5*R, -(25*H/2));

  rho1 = RHO1;
  rho2 = rho1 / RHOR;
  f.sigma = SIGMA;

  const scalar c[] = theta0*pi/180.;
  contact_angle = c;

  TOLERANCE = 1e-4 [*];

  mu1 = MUl;
  mu2 = MUg;
  run(); 

}

event init (t = 0) {

  // Refine interfaces until they reach the maximum level
  scalar ibmFake[], fFake[];
  face vector ibmfFake[];
  astats ss;
  int ic = 0;
  do {
    ic++;
    lq_domain (fFake);
    solid_domain (ibmFake, ibmfFake);
    ss = adapt_wavelet ({ibmFake, fFake}, (double[]) {1.e-30, 1.e-30},
		        maxlevel = LEVEL , minlevel = (1)); 
  } while ((ss.nf || ss.nc) && ic < 1000);
  
  // Reinitalize the solid and fluid/gas fields
  lq_domain (f);
  solid_domain (ibm, ibmf);

  // Initialize velocity
  foreach() {
    foreach_dimension()
      u.x[] = 0.;
  }
  boundary({u, f});
}

scalar pid[];

event logfile (i++; t <= t_end) {

  fprintf (stderr, "%d %g\n", i, t);

}


event snapshot (t == t_end) {
  scalar omega[];
  vorticity (u, omega);
  foreach()
    omega[] *= H/Ug;

  char name[80];
  sprintf (name, "%d-vort-%g.png", k, t);
  FILE * fp1 = fopen (name, "w");
  view (fov = 7, width = 1200, height = 600);
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f", filled = 1, fc = {1,0,0}, lw = 5);
  squares ("omega", max = 10, min = -10, map = cool_warm);
  save (fp = fp1);
  fclose (fp1);
  
  scalar umag[];
  foreach()
    umag[] = f[] * sqrt(sq(u.x[]) + sq(u.y[])) / Ug;

  sprintf (name, "%d-%d-velo-%g.png", k, j, t);
  FILE * fp3 = fopen (name, "w");
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f");
  squares ("umag", min = 0);
  save (fp = fp3);
  fclose (fp3);

  sprintf (name, "%d-%d-grid-%g.png", k, j, t);
  FILE * fp4 = fopen (name, "w");
  clear();
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f");
  cells();
  save (fp = fp4);
  fclose (fp4);
}

event movie (t += 0.05; t <= t_end) {

  char name[80];
  sprintf (name, "%d-%d-vof-%g.png",k, j, t);
  FILE * fp1 = fopen (name, "w");
  view (fov = 7, width = 1200, height = 600);
  draw_vof ("ibm", "ibmf", filled = -1);
  draw_vof ("f");
  squares ("f", min = 0, max = 1);
  save (fp = fp1);
  fclose (fp1);
}

#if 1
event adapt (i++) {
  scalar ibm1[], f1[];
  foreach() {
    ibm1[] = ibm[];
    f1[] = f[];
  }
  adapt_wavelet ({ibm1,f,u}, (double[]){1e-1,1e-1,uemax,uemax},
		  maxlevel = LEVEL, minlevel = LEVEL - 4);
}
#endif

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", LEVEL, s.real, i, s.speed);
  fprintf (fp, "%d\t%d\t%d\t%d\n", mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  fflush (fp);
  return 1;
}
