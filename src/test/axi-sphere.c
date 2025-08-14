#include "../ibm-gcm.h"
//#include "../my-axi.h"
#include "../my-centered.h"
//#include "navier-stokes/centered.h"
#include "../ibm-gcm-events.h"
#include "navier-stokes/perfs.h"

#define L0 2.
#define D 0.5
#define LEVEL 9
#define U0 1.

#define RE 100

u.n[left] = dirichlet(U0);
p[left] = neumann(0);
pf[left] = neumann(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

face vector muv[];

int main() {

  size (L0);
  init_grid (1 << LEVEL);
  origin(-L0/2., 0);
 
  //stokes = true;

  mu = muv;

  DT = 2e-3;
  run();
}

event init (t = 0) {

  solid (ibm, ibmf, sq(x) + sq(y) - sq(D/2.));
  event("update_metric");

  foreach()
    u.x[] = ibm[]*U0;

}

face vector alphav[];
scalar rhov[], pid[];
#if 1
event properties (i++) {
  foreach_face() {
    muv.x[] = fm.x[]*(U0*D)/RE;
    alphav.x[] = alpha.x[];
  }
  foreach() {
    rhov[] = rho[];
    pid[] = pid();
  }
  boundary({muv});
}
#endif

event logfile (i ++; t <= 5) {

    foreach_face() 
      if (!fm.x[])
        uf.x[] = 0;
  boundary({uf});

  coord Fp = {0}, Fmu = {0};
  ibm_force(p, u, muv, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x) / (0.5 * (U0) * sq(D/2.) * M_PI);
  double CL = (Fp.y + Fmu.y) / (0.5 * (U0) * sq(D/2.) * M_PI);

  fprintf (stderr, "%d %g %g %g %g %g\n", i, t, 
        normf(u.x).max, normf(u.y).max, CD, CL);
}

#if TREE && 1
event adapt (i++) {
  adapt_wavelet ({ibm, u.x, u.y}, (double[]){1e-1, 1e-1, 1e-1},
		 maxlevel = LEVEL, minlevel = LEVEL - 5);
  event("update_metric");
}
#endif

