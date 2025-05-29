#include "fractions.h"


attribute {
  scalar * tracers, c;
  bool inverse;
}


extern scalar * interfaces;
extern face vector uf;
extern double dt;

double cerror0, cerror1_x, cerror1_y, cerror2;

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
scalar c0[];

vector nfg[], nsg[];
scalar alphafg[], alphasg[];

scalar flux[];  // real flux
scalar fluxf[]; // fake flux

foreach_dimension()
static void sweep_x (scalar c, scalar cc, scalar * tcl, scalar cr, scalar ccr, scalar ibm0, face vector ibmf0)
{
  vector nf[], ns[];
  scalar alphaf[], alphas[];
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
  reconstruction (c, nf, alphaf);
  reconstruction (ibm0, ns, alphas);
  //immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);
  //reconstruction (c, nf, alphaf);
  //reconstruction (ibm0, ns, alphas);
  cerror1_x = real_volume (c);

  foreach_face(x, reduction (max:cfl)) {

#if IBM
    double un = uf.x[]*dt/(Delta), s = sign(un);
#else
    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
#endif
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
    if (cs[] >= 1.)
#elif IBM
    if (ibm0[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    double cf = 0;
    coord tempnf = {-s*nf.x[i], nf.y[i], nf.z[i]};
    coord lhs = {-0.5, -0.5, -0.5}, rhs = {s*un - 0.5, 0.5, 0.5};

    if (ibm0[i] >= 1. || ibm0[i] <= 0.) {
        cf = (c[i] <= 0. || c[i] >= 1.)? c[i] : rectangle_fraction (tempnf, alphaf[i], lhs, rhs);
    }
    else if (ibm0[i] > 0. && ibm0[i] < 1.) {
        coord tempns = {-s*ns.x[i], ns.y[i], ns.z[i]};
        if (cr[i] <= 0. || un <= 0.)
            cf = 0.;
        else if (cr[i] >= ibm0[i]) { // interfacial cell is full
            double alphac = plane_alpha (cr[i], (coord){ns.x[i], ns.y[i]});
            cf = (cr[i] <= 0. || cr[i] >= 1.)? cr[i] : rectangle_fraction (tempns, alphac, lhs, rhs);
        }
        else {
            cf = immersed_fraction (c[i], tempnf, alphaf[i], tempns, alphas[i], lhs, rhs,0);
            #if 0
            fprintf(stderr, "VOF (b): %g (%g, %g) ibm[%d]=%g cf=%0.15g"
                            " c[%d]=%0.15g cr[%d]=%0.15g un=%g, uf=%g\n", 
                             indicator.x, x, y, i, ibm[i], cf, i, c[i], i, cr[i], 
                             un, uf.x[]);
           #endif
       }
    }
    else {
        cf = rectangle_fraction (tempnf, alphaf[i], lhs, rhs);
    }

    flux[] = cf*ibmf0.x[]*uf.x[];

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
 
  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "src/vof.h:%d: warning: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     __LINE__, cfl - 0.5), fflush (ferr);

  foreach()
    if (ibm0[] > 0) {

#if 1
      if (on_interface(ibm0)) 
        fprintf (stderr, "F* W/O FLUX: %g (%g, %g) c[]=%0.15g cr[]=%0.15g"
                         " ibmf[]=%g ibmf[1]=%g uf[]=%g uf[1]=%g dx=%g\n",
                         indicator.x, x, y, c[], cr[], ibmf.x[], ibmf.x[1],
                         uf.x[], uf.x[1], Delta);
#endif
      if (on_interface(ibm0)) {
          coord mp, n;
          //double area = ibm0_geometry (point, &mp, &n);
          double area = ibm_geometry (point, &mp, &n);
          double mpx, mpy, mpz;
          local_to_global (point, mp, &mpx, &mpy, &mpz);

          double divs = uibm_x(mpx,mpy,mpz) * n.x * area;

          c[]  += dt*(flux[] - flux[1] + cc[]*(ibmf0.x[1]*uf.x[1] - ibmf0.x[]*uf.x[] - divs))/Delta;
          cr[] += dt*(flux[] - flux[1] + ccr[]*(ibmf0.x[1]*uf.x[1] - ibmf0.x[]*uf.x[] - divs))/(ibm[]*Delta);
      }
      else {
          c[]  += dt*(flux[] - flux[1] + cc[]*(ibmf0.x[1]*uf.x[1] - ibmf0.x[]*uf.x[]))/(Delta);
          cr[] += dt*(flux[] - flux[1] + ccr[]*(ibmf0.x[1]*uf.x[1] - ibmf0.x[]*uf.x[]))/(Delta);
      }
#if 1
      if (on_interface(ibm)) 
        fprintf (stderr, "F* W/FLUX: %g (%g, %g) c[]=%0.15g cr[]=%0.15g cc[]=%g"
                         " flux[]=%g flux[1]=%g div=%g\n",
                         indicator.x, x, y, c[], cr[], cc[], flux[], flux[1], divg1[]);
#endif
    }

#if 0
  reconstruction (c, nf, alphaf);
  reconstruction (ibm, ns, alphas);
  immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);
#endif

  delete (tfluxl); free (tfluxl);
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {
    vector nf[], ns[];
    scalar alphaf[], alphas[], creal[];

    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */

    scalar cc[], ccr[], * tcl = NULL, * tracers = c.tracers;
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
      ccr[] = (c[] > 0.5*ibm[]);
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

   reconstruction (c, nf, alphaf);
   reconstruction (ibm, ns, alphas);

   // TODO: cr should be = 1 in cells where cr == ibm, i.e. full cells,
   //       so i think we just dont multiply by ibm[] here.
   foreach() {
     if (on_interface(ibm) && on_interface(c)) {
       creal[] = immersed_fraction (c[], (coord){nf.x[], nf.y[]}, alphaf[],
                                         (coord){ns.x[], ns.y[]}, alphas[],
                                         (coord){-0.5, -0.5, -0.5},
                                         (coord){0.5, 0.5, 0.5}, 0);
     }
     else
       creal[] = c[];
     
     cr0[] = creal[];   
     c0[] = c[];
   }

    cerror0 = real_volume (c);

    void (* sweep[dimension]) (scalar, scalar, scalar *, scalar, scalar, scalar, face vector);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++) {
      char ind = (i + d) % dimension == 0? 'x': 'y';
      fprintf(stderr, "\n=== %c SWEEP (i = %d) ===\n", ind, i);
      sweep[(i + d) % dimension] (c, cc, tcl, creal, ccr, ibm, ibmf);
    }
    delete (tcl), free (tcl);

   foreach() {
     cr[] = creal[];
   }

   cerror2 = real_volume (c);

   reconstruction (c, nf, alphaf);
   reconstruction (ibm, ns, alphas);
   immersed_reconstruction (c, creal, nf, alphaf, ns, alphas);

   reconstruction (c, nfg, alphafg);
   reconstruction (ibm, nsg, alphasg);
  }
}

event vof (i++)
  vof_advection (interfaces, i);

