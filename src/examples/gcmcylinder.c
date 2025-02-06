#include "../ibm-gcm.h" 
#include "../my-centered.h"
#include "../ibm-gcm-events.h" 
#include "view.h"

#define L0 15.
#define D 0.5
#define LEVEL 10

int Re;
int minlevel = LEVEL - 5;
double U0 =  1.;             // inlet velocity
double t_end = 50;
double tf_start = 25;
coord ci = {L0/4., L0/2.};     // initial coordinates of cylinder

face vector muv[];


u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

u_ibm_dirichlet_x(x*x + y*y > 1? 0: 0)
u_ibm_dirichlet_y(x*x + y*y > 1? 0: 0)

int main() {
  size(L0);
  init_grid (1 << (LEVEL - 2));
  mu = muv;
  TOLERANCE = 1.e-6;
  CFL = 0.8;
  
  Re = 100;
  run();

}

event init (t = 0) {
  //initial_refine (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) - sq(D/2), LEVEL,minlevel,1000); 
  solid (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) - sq(D/2));
  refine (ibm[] < 1 && level < LEVEL);
  solid (ibm, ibmf, sq(x - ci.x) + sq(y - ci.y) - sq(D/2));

  // Must initialize velocity field to prevent noisy start
  foreach() {
    u.x[] = ibm[] * U0;
    u.y[] = 0;
  }

}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

event logfile (i++) {
  coord Fp, Fmu;
  ibm_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x) / (0.5*sq(U0)*(D));
  double CL = (Fp.y) / (0.5*sq(U0)*(D));

  double CD1 = (Fmu.x) / (0.5*sq(U0)*(D));
  double CL1 = (Fmu.y) / (0.5*sq(U0)*(D));


  fprintf (stderr, "%d %g %d %d %d %d %d %d %d %g %g %g %g\n",
          i, t, Re, mgpf.i, mgpf.nrelax, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax,
          CD, CL, CD1, CL1);
}

event rings (t = end)
{
  char name[80];
  scalar ring[];
  fraction (ring, - sq(x - ci.x) - sq(y - ci.y ) + sq((0.5+0.0146)/2));
  sprintf (name, "ring1-%d", Re);
  FILE * fv6 = fopen (name, "w");
  foreach()
    if (ring[] > 0 && ring[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv6, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv6);
  fclose (fv6);

  fraction (ring, - sq(x - ci.x ) - sq(y - ci.y ) + sq((0.5+0.0293)/2));
  sprintf (name, "ring2-%d", Re);
  FILE * fv7 = fopen (name, "w");
  foreach()
    if (ring[] > 0 && ring[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv7, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv7);
  fclose (fv7);

  fraction (ring, - sq(x - ci.x ) - sq(y - ci.y ) + sq((0.5+0.0439)/2));
  sprintf (name, "ring3-%d", Re);
  FILE * fv8 = fopen (name, "w");
  foreach()
    if (ring[] > 0 && ring[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv8, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv8);
  fclose (fv8);
 
}

int count = 0;
event frequency (i++) {
  if (t >= tf_start && Re >= 50) {
    char name[80];
    sprintf (name, "freq-7.5-%d.dat", Re);
    FILE * fp = fopen (name, "a");
    foreach_point(7.5000, 5.0000) {
      fprintf (fp, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
    }
    fclose (fp);

    sprintf (name, "freq-5-%d.dat", Re);
    FILE * fp1 = fopen (name, "a");
    foreach_point(5.0000, 5.0000) {
      fprintf (fp1, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
    }
    fclose (fp1);
    count++;
  }
}


event profile (t = t_end) {
  double delta = L0/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx2-%d", Re); // x = center of cylinder
  FILE * fv1 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (ci.x, i) {
      fprintf (fv1, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv1);
  fclose (fv1);

  sprintf (name, "vprofx3-%d", Re); // x = outlet
  FILE * fv2 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (L0-delta, i) {
      fprintf (fv2, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]); 
    }
  }
  fflush (fv2);
  fclose (fv2);

  sprintf (name, "vprofy1-%d", Re); // y = L0/2
  FILE * fv3 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (i, ci.y) {
      fprintf (fv3, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv3);
  fclose (fv3);
}

event adapt (i++) {
  scalar ibmsf[]; 
  foreach() {
    ibmsf[] = vertex_average (point, ibm);
  }
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,3e-3,3e-3},
                           maxlevel = LEVEL, minlevel = minlevel);
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fprintf (fp, "%d\t%d\t%d\t%d\n", mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  fflush (fp);
  return 1;
}

