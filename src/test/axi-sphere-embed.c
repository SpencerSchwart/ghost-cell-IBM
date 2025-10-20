#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

#define L0 2.
#define D 0.5
#define LEVEL 8
#define U0 1.

int RE = 250;

double tend = 10;

u.n[left] = dirichlet(U0);
p[left] = neumann(0);
pf[left] = neumann(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);

u.t[embed] = dirichlet(0);
u.n[embed] = dirichlet(0);

face vector muv[];

int main() {

  size (L0);
  init_grid (1 << LEVEL);
  origin(-L0/2., 0);
 
  //stokes = true;

  TOLERANCE=1e-6;

  mu = muv;

  run();
}

event init (t = 0) {

  solid (cs, fs, sq(x) + sq(y) - sq(D/2.));

  restriction({cs});

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);

  cm.refine = cm.prolongation = refine_cm_axi;
  cs.refine = cs.prolongation = fraction_refine;
  fm.x.refine = refine_face_x_axi;
  fm.y.refine = refine_face_y_axi;
  metric_embed_factor = axi_factor;

  restriction ({cs, fs, cm, fm});

  foreach()
    u.x[] = cs[]*U0;

}

event properties (i++) {
  foreach_face() {
    muv.x[] = fm.x[]*(U0*D)/RE;
  }
  boundary({muv});
}

event logfile (i ++; t <= tend) {

    foreach_face() 
      if (!fm.x[])
        uf.x[] = 0;
  boundary({uf});

  coord Fp = {0}, Fmu = {0};
  embed_force(p, u, muv, &Fp, &Fmu);

  double CD = 2*pi*(Fp.x + Fmu.x) / (0.5 * (U0) * sq(D/2.) * M_PI);
  double CL = 2*pi*(Fp.y + Fmu.y) / (0.5 * (U0) * sq(D/2.) * M_PI);

  fprintf (stderr, "%d %g %g %g %g %g\n", i, t, 
        normf(u.x).max, normf(u.y).max, CD, CL);
  fprintf (stdout, "%d %g %g %g %g %g\n", i, t, 
        normf(u.x).max, normf(u.y).max, CD, CL);
}

#if 0
event surf_data (t = end)
{
  char name[80];
  sprintf (name, "%d-surface", RE);
  FILE * fpxy = fopen (name, "w");

  scalar cf[];
  //skin_friction (u, mu, cf);

  foreach () {
    if (cs[] > 0 && cs[] < 1) {
      coord midPoint, n, b, cc = {x,y,z};
      double area = embed_geometry (point, &b, &n);

      foreach_dimension()
        midPoint.x = cc.x + b.x*Delta;

      double theta = atan2(y, x) * (180/M_PI) + 180;
      double cp = extrapolate_scalar (point, ibm, midPoint, n, p) / (0.5 * sq(U0));
      double cfi = cf[] / (0.5 * sq(U0));

      fprintf (fpxy, "%g %g %g %g %g %g\n", 
                      x, y, theta, area, cp, cfi);
    }
  }
  fclose(fpxy);
}
#endif

#if 1
event adapt (i++) {
  adapt_wavelet ({cs, u.x, u.y}, (double[]){1e-3, 1e-3, 1e-3},
		 maxlevel = LEVEL, minlevel = LEVEL - 5);

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
}
#endif
