/**
# Impact of a viscoelastic drop on a solid

We solve the axisymmetric, incompressible, variable-density,
Navier--Stokes equations with two phases and use the log-conformation
method to include viscoelastic stresses. The curvature module is used
to compute interface properties. */

#include "../ibm-gcm.h"
//#include "embed.h"
//#include "navier-stokes/centered.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../my-two-phase.h"
//#include "tension.h"
#include "curvature.h"

/**
The density and viscosity ratios are 1000. The Reynolds number based
on the droplet diameter, velocity and viscosity is 5 and the Froude
number is 2.26. */

#define RHO_r 0.001
#define MU_r 0.001
#define RE 5.
#define FR 2.26
#define LEVEL 7

/**
The dimensionless viscoelastic properties used are the ratio of
the solvent to the total viscoelastic viscosity (polymeric
plus solvent) and the Weissenberg number. */

#define BETA 0.1
#define WI 1.0

/**
The drop comes from the right. We allow the fluid to get through that
boundary. */

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

/**
The wall is at the left side. We apply a no-slip boundary condition and a 
non-wetting condition for the VOF tracer. */

u.t[left] = dirichlet(0);
f[left]   = 0.;

#if IBM
u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (0)
#endif

int main() {

  /**
  The domain spans $[0:2.6]\times[0:2.6]$. */

  size (2.6);
  init_grid (1 << LEVEL);
  
  /**
  The densities and viscosities are defined by the parameters above. */

  rho1 = 1.;
  rho2 = RHO_r;
  mu1 = BETA/RE;
  mu2 = MU_r/RE;
  // f.sigma = 1;

  /**
  The viscoelastic fields will be set below. */
  
  /**
  We set a maximum timestep. This is necessary for proper temporal
  resolution of the viscoelastic stresses. */
  
  DT = 2e-3;
  run();
}

event init (t = 0) {

  /**
  The drop is centered on (2,0) and has a radius of 0.5. */

  fraction (f, - sq(x - 2.) - sq(y) + sq(0.5));

#if IBM
  solid (ibm, ibmf, x - 0.9825);
#endif
#if EMBED
  solid (cs, fs, x - 0.9825);
  u.t[embed] = dirichlet(0);
  u.n[embed] = dirichlet(0);
  f[embed] = dirichlet(0);
#endif
  /**
  The initial velocity of the droplet is -1. */

  foreach()
    u.x[] = - f[];
}


/**
We add the acceleration of gravity. */
#if 0
event acceleration (i++) {
  face vector av[];
  foreach_face(x)
    av.x[] -= 1./sq(FR);
  av = a;
}
#endif

/**
We adapt the solution at every timestep based on the interface and
velocity errors. */

#if 0
event adapt (i++) {
  adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 5e-3, 5e-3},
		 maxlevel = LEVEL, minlevel = LEVEL - 2);
}
#endif

/**
We track the spreading diameter of the droplet. */

event logfile (i += 20; t <= 5) {
  scalar pos[];
  position (f, pos, {0,1});
  fprintf (stderr, "%g %g\n", t, 2.*statsf(pos).max);
}

/**
We generate a movie of the interface shape. */
#if 0
#include "view.h"

event viewing (i += 10) {
  view (width = 400, height = 400, fov = 20, ty = -0.5,
	quat = {0, 0, -0.707, 0.707});

  clear();
  draw_vof ("f", lw = 2);
  squares ("u.x", linear = true);
  box (notics = true);
  mirror ({0,1}) {
    draw_vof ("f", lw = 2);
    squares ("u.y", linear = true);
    box (notics = true);
  }
  save ("movie.mp4");
  
  /**
  We can optionally visualise the results while we run. */

#if 0
  static FILE * fp = popen ("bppm","w");
  save (fp = fp);
#endif
}
#endif

