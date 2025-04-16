double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0)
{
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  
  if (mu1 || mu2)
    mu = new face vector;

  /**
  We add the interface to the default display. */

  display ("draw_vof (c = 'f');");
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#if FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++)
{
  
  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif // !sf

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

#include "fractions.h"

event properties (i++)
{

#if IBM
#if 0
  vector nf[], ns[];
  scalar alphaf[], alphas[];
  reconstruction (f, nf, alphaf);
  reconstruction (ibm, ns, alphas);
  
  foreach_face() {
    coord nfluid0, nsolid0, nfluid1, nsolid1;
      nfluid0.x = nf.x[];
      nsolid0.x = ns.x[];
      nfluid1.x = nf.x[-1];
      nsolid1.x = ns.x[-1];
      nfluid0.y = nf.y[];
      nsolid0.y = ns.y[];
      nfluid1.y = nf.y[-1];
      nsolid1.y = ns.y[-1];

    double freal0 = sf[], freal1 = sf[-1];
    //double ff = (sf[] + sf[-1])/2.;
    if (on_interface(ibm) && on_interface(f)) {
        freal0 = immersed_fraction (f[], nfluid0, alphaf[], nsolid0, alphas[],
                                           (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);
    }
    else if (on_interface(ibm) && f[] >= 1.-1e-6) {
        freal0 = f[]*ibm[];
    }
    if (ibm[-1] > 0 && ibm[-1] < 1 && f[-1] < 1 && f[-1] > 0) {
        freal1 = immersed_fraction (f[-1], nfluid1, alphaf[-1], nsolid1, alphas[-1],
                                           (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);
    }
    else if (ibm[-1] > 0 && ibm[-1] < 1 && f[-1] >= 1.-1e-6) {
        freal1 = f[-1]*ibm[-1];
    }
    double ff = (freal0 + freal1)/2.; 
    alphav.x[] = ibmf.x[]/(rho(ff)+SEPS);
    if (mu1 || mu2) {
      face vector muv = mu;
      muv.x[] = fm.x[]*mu(ff); // should fm be ibmf here?
    }
  }
#endif
    foreach_face() {
        double ff = (sf[] + sf[-1])/2.;
        alphav.x[] = ibmf.x[]/rho(ff);
        if (mu1 || mu2) {
          face vector muv = mu;
           muv.x[] = fm.x[]*mu(ff); // should fm be ibmf here?
        }
    }
#else // !IBM
    foreach_face() {
        double ff = (sf[] + sf[-1])/2.;
        alphav.x[] = fm.x[]/rho(ff);
        if (mu1 || mu2) {
          face vector muv = mu;
           muv.x[] = fm.x[]*mu(ff); // should fm be ibmf here?
        }
    }
#endif // !IBM
  
  foreach() {
#if IBM
    rhov[] = ibm[]*rho(sf[]);
#else
    rhov[] = cm[]*rho(sf[]);
#endif
  }

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}
