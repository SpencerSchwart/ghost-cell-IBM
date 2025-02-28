#include "fractions.h"


attribute {
  scalar * tracers, c;
  bool inverse;
}


extern scalar * interfaces;
extern face vector uf;
extern double dt;


foreach_dimension()
static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
	if (t.gradient)
	  return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
	else
	  return (t[1]/cr - t[-1]/cl)/(2.*Delta);
      }
      else
	return (t[1]/cr - t[]/cc)/Delta;
    }
    else if (cl >= cmin)
      return (t[]/cc - t[-1]/cl)/Delta;
  }
  return 0.;
}


#if TREE
static void vof_concentration_refine (Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = 0.;
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
    double sc = s.inverse ? s[]/(1. - f[]) : s[]/f[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
	s[] += child.x*g.x*cm[-child.x]/cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}


event defaults (i = 0)
{
  for (scalar c in interfaces) {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
    }
  }
}
#endif // TREE


event defaults (i = 0)
{
  for (scalar c in interfaces) {
    scalar * tracers = c.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, c);
  }
}


event stability (i++) {
  if (CFL > 0.5)
    CFL = 0.5;
}

coord indicator = {0,1};

scalar cr[]; // holds the real fluid
scalar cr0[];

vector nf[], ns[];
scalar alphaf[], alphas[];

void reconstruction_contact_vof (scalar f, vector n, scalar alpha);

foreach_dimension()
static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  //vector n[], ns[];
  //scalar alpha[], alphas[], flux[];
  scalar flux[];  // real flux
  scalar fluxf[]; // fake flux
  double cfl = 0.;


  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }


    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl)
	gf[] = vof_concentration_gradient_x (point, c, t);
    }
  }
  
  // 1. Find n and alpha for f[] and ibm[]
  
  trash({ns,nf,alphaf,alphaf});

  reconstruction (c, nf, alphaf);
  //reconstruction_contact_vof (c, nf, alphaf);
  reconstruction (ibm, ns, alphas);

  foreach_face(x, reduction (max:cfl)) {

#if IBM
    double un = uf.x[]*dt/(Delta + SEPS), s = sign(un);
#else
    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
#endif
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
    if (cs[] >= 1.)
#elif IBM
    if (ibm[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

#if 0
    double cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
      rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			  (coord){-0.5, -0.5, -0.5},
			  (coord){s*un - 0.5, 0.5, 0.5});
#endif
#if 1
    double cf = 0;
    //if ((c[i] <= 0. || c[i] >= 1.) && (ibm[i] <= 0. || ibm[i] >= 1))
    if ((c[i] <= 0. || c[i] >= 1.)) {
        cf = c[i];
        fluxf[] = 0;
    }
#if 1
    else if (ibm[i] > 0. && ibm[i] < 1. && c[i] > 0) {
        coord tempnf = {-s*nf.x[i], nf.y[i], nf.z[i]};
        coord tempns = {-s*ns.x[i], ns.y[i], ns.z[i]};
        coord lhs = {-0.5, -0.5, -0.5}, rhs = {s*un - 0.5, 0.5, 0.5};
        cf = immersed_fraction (c[i], tempnf, alphaf[i], tempns, alphas[i], lhs, rhs,0);
        if (cf == 0)
            fluxf[] = uf.x[]*ibmf.x[]*rectangle_fraction (tempnf, alphaf[i], lhs, rhs);
        else
            fluxf[] = 0;
    }
#endif
    else {
        coord tempnf = {-s*nf.x[i], nf.y[i], nf.z[i]};
        coord lhs = {-0.5, -0.5, -0.5}, rhs = {s*un - 0.5, 0.5, 0.5};
        cf = rectangle_fraction (tempnf, alphaf[i], lhs, rhs);
        fluxf[] = 0;
    }
    /*
    double cf = ((c[i] <= 0. || c[i] >= 1.) && (ibm[i] <= 0. || ibm[i] >= 1)) ? 
                c[i] : (ibm[i] <= 0. || ibm[] >= 1.)?
      rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			  (coord){-0.5, -0.5, -0.5},
			  (coord){s*un - 0.5, 0.5, 0.5}):
      immersed_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
                         (coord){-s*ns.x[i], ns.y[i], ns.z[i]}, alphas[i],
                         (coord){-0.5, -0.5, -0.5},
                         (coord){s*un - 0.5, 0.5, 0.5});
    */
#endif
#if 0
    if (ibm[i] > 0 && ibm[i] < 1) {
      fprintf(stderr, "\n### %g ibm[i] = %g f[i]=%g cf=%g ###\n", indicator.x, ibm[i], c[i], cf);
        double f0 = immersed_fraction ((coord){-s*nf.x[i], nf.y[i], nf.z[i]}, alphaf[i],
                         (coord){-s*ns.x[i], ns.y[i], ns.z[i]}, alphas[i],
                         (coord){-0.5, -0.5, -0.5},
                         (coord){s*un - 0.5, 0.5, 0.5});
        double alpha0 = immersed_line_alpha ((coord){-s*nf.x[i], nf.y[i], nf.z[i]}, alphaf[i],
                         (coord){-s*ns.x[i], ns.y[i], ns.z[i]}, alphas[i],f0);
        fprintf(stderr, "alpha0 = %g vs alphaf = %g | f0=%g cf=%g\n", alpha0, alphaf[i], f0, cf);
    }
#endif


    /**
    Once we have the upwind volume fraction *cf*, the volume fraction
    flux through the face is simply: */

    flux[] = cf*ibmf.x[]*uf.x[]; // add ibmf.x[] to uf.x[]?
    //flux[] = cf*uf.x[]; // add ibmf.x[] to uf.x[]?

    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
    	cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	    double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
    	tflux[] = ff*cf1*uf.x[];
      }
      else
	    tflux[] = 0.;
    }
  }
  delete (gfl); free (gfl);
 
  foreach() {
    if (ibm[] > 0 && ibm[] < 1 && c[] > 0 && c[] < 1) {
      //fprintf(stderr, "\n### %g ibm[i] = %g f[i]=%g cf=%g ###\n", indicator.x, ibm[i], c[i], cf);
        // get the real volume fraction of interfacial cells before advection
        cr[] = immersed_fraction (c[], (coord){nf.x[], nf.y[]}, alphaf[],
                         (coord){ns.x[], ns.y[]}, alphas[],
                         (coord){-0.5, -0.5, -0.5},
                         (coord){0.5, 0.5, 0.5}, 0);
        //fprintf(stderr, "cr0[] = %g f[] = %g ibm[]=%g\n", cr[], c[], ibm[]);
    }
    else {
        cr[] = ibm[]*c[];
        //fprintf(stderr, "cr0[] = %g f[] = %g ibm[]=%g\n", cr[], c[], ibm[]);
    }
    cr0[] = cr[];
  }
  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "src/vof.h:%d: warning: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     __LINE__, cfl - 0.5), fflush (ferr);



  foreach()
    if (ibm[] > 0.) {
      c[] += dt*((flux[] + fluxf[]) - (flux[1] - fluxf[1]) + cc[]*(ibmf.x[1]*uf.x[1] - ibmf.x[]*uf.x[]))/(Delta);
      //c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/(Delta);
      // add the real flux to the inital real volume
      cr[] += dt*(flux[] - flux[1] + cc[]*(ibmf.x[1]*uf.x[1] - ibmf.x[]*uf.x[]))/(Delta);

      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
    	t[] += dt*ibm[]*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(ibm[]*Delta);
    }

  immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);

  delete (tfluxl); free (tfluxl);
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {

    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */

    scalar cc[], * tcl = NULL, * tracers = c.tracers;    
    for (scalar t in tracers) {
#if !NO_1D_COMPRESSION
      scalar tc = new scalar;
      tcl = list_append (tcl, tc);
#endif // !NO_1D_COMPRESSION
#if TREE
      if (t.refine != vof_concentration_refine) {
	t.refine = t.prolongation = vof_concentration_refine;
	t.restriction = restriction_volume_average;
	t.dirty = true;
	t.c = c;
      }
#endif // TREE
    }
    foreach() {
      cc[] = (c[] > 0.5);
#if !NO_1D_COMPRESSION
      scalar t, tc;
      for (t, tc in tracers, tcl) {
	if (t.inverse)
	  tc[] = c[] < 0.5 ? t[]/(1. - c[]) : 0.;
	else
	  tc[] = c[] > 0.5 ? t[]/c[] : 0.;
      }
#endif // !NO_1D_COMPRESSION
    }

    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */

    void (* sweep[dimension]) (scalar, scalar, scalar *);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++)
      sweep[(i + d) % dimension] (c, cc, tcl);
    delete (tcl), free (tcl);

  //reconstruction (c, nf, alphaf);
  //reconstruction (ibm, ns, alphas);
  //immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);
  }
}

event vof (i++)
  vof_advection (interfaces, i);

/**
## References

~~~bib
@Article{lopez2015,
  title = {A VOF numerical study on the electrokinetic effects in the 
           breakup of electrified jets},
  author = {J. M. Lopez-Herrera and A. M. Ganan-Calvo and S. Popinet and
            M. A. Herrada},
  journal = {International Journal of Multiphase Flows},
  pages = {14-22},
  volume = {71},
  year = {2015},
  doi = {doi.org/10.1016/j.ijmultiphaseflow.2014.12.005},
  url = {http://gerris.dalembert.upmc.fr/papers/lopez2015.pdf}
}
~~~
*/
