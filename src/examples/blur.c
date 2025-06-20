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

double MUl[4] = {2.e-2, 5.e-3, 5.e-3, 5.e-3};
double MUg[4] = {1.e-1, 2.5e-2, 2.5e-2, 2.5e-2};
double Ul[4] = {20, 2, 1, 0.5};
double Ug[4] = {40, 160, 160, 160};
double theta0[3] = {70, 30, 90};

u.n[left] = dirichlet (fabs(y) < R? Ul[j]:0);
u.t[left] = dirichlet (0);
p[left] = neumann (0);
pf[left] = neumann (0);

u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

u.n[top] = dirichlet (x > 0 && x <= H? -Ug[j]:0);
u.t[top] = dirichlet (0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = dirichlet (x > 0 && x <= H? Ug[j]:0);
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
    // if (x < -H && fabs(y) < R)
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

  const scalar c[] = theta0[0]*pi/180.;
  contact_angle = c;

  TOLERANCE = 1e-4 [*];
#if 0
  for (j = 0; j <= 3; j++) {
    mu1 = MUl[j];
    mu2 = MUg[j];
    // DT = 0.1;
    run();
  }
#endif
    mu1 = MUl[3];
    mu2 = MUg[3];
    run(); 
  /*
  for (k = 1; k <= 2; k++) {
    const scalar cc[] = theta0[k]*pi/180;
    contact_angle = cc;
    run();
  }
  */
}

event init (t = 0) {
  // mask(y > 5*H ? top: y < -5*H ? bottom : none);

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
  //refine (ibm[] > 0 && ibm[] < 1 && level < LEVEL);
  
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
  double minFM = 1e30, maxUF=-1e-30;
  foreach_face(reduction(min:minFM) reduction(max:maxUF)) {
    if (fm.x[] && fm.x[] < minFM)
        minFM = fm.x[];
    if (fabs(uf.x[]) > fabs(maxUF))
        maxUF = uf.x[];
  }

  foreach()
    pid[] = pid();

  fprintf (stderr, "%d %g %d %d %g %g\n", i, t, j, k, minFM, maxUF);

  /*
  char name[80];
  sprintf (name, "%d-%d-xh2", k, j);
  FILE * fp2 = fopen (name, "a");
  foreach_point(-2*H, 0)
    fprintf (fp2, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp2);

  sprintf (name, "%d-%d-xh1", k, j);
  FILE * fp1 = fopen (name, "a");
  foreach_point(-1*H, 0)
    fprintf (fp1, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp1);

  sprintf (name, "%d-%d-xh0", k, j);
  FILE * fp0 = fopen (name, "a");
  foreach_point(0, 0)
    fprintf (fp0, "%g %g %g %g %g\n", t, x, y, f[], u.y[]);
  fclose(fp0);
  */
}

/*
event interface (t += 0.3125; t <= t_end) {
  char name[80];
  sprintf (name, "%d-%d-int-%g", k, j, t);
  FILE * fpi = fopen(name, "w");
  output_facets (f, fp = fpi);
  fclose (fpi);
}
*/

event snapshot (t == t_end) {
  scalar omega[];
  vorticity (u, omega);
  foreach()
    omega[] *= H/Ug[j];

  char name[80];
  sprintf (name, "%d-%d-vort-%g.png", k, j, t);
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
    umag[] = f[] * sqrt(sq(u.x[]) + sq(u.y[])) / Ug[j];

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
    //ibm1[] = 0;
    f1[] = f[];
  }
  adapt_wavelet ({ibm,ibm1,f,u}, (double[]){1e-1,1e-1,1e-1,uemax,uemax},
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
