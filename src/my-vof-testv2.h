#include "fractions.h"

#define VPRINT 0

// TODO: use ibmf to calculate ibm normal during reconstruction!!!
//      - reconstruction uses MYC, while ibm_geomtry uses ibmf to find n

attribute {
  scalar * tracers, c;
  bool inverse;
}

scalar cid[];

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

scalar ch[], cr[];
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
  ch.refine = ch.prolongation = fraction_refine;
  ch.dirty = true;
  cr.refine = cr.prolongation = fraction_refine;
  cr.dirty = true;
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

void set_contact_angle (scalar f, const scalar fr0, const scalar ibm, 
                        vector nu, scalar alphau, vector ns, scalar alphas, scalar id);
void set_contact_angle_tension (scalar f, const scalar fr0, const scalar ibm,
                                vector nf, scalar alphaf, vector ns, scalar alphas);
coord indicator = {0,1,2};

vector nfg[], nsg[], divs[];
scalar alphafg[], alphasg[];

face vector ibmf_temp[];
vector fluxr[];

void move_solid_x(scalar ibm, face vector ibmf);
void move_solid_y(scalar ibm, face vector ibmf);
void move_solid_z(scalar ibm, face vector ibmf);

foreach_dimension()
static void sweep_x (scalar c, scalar ch, scalar cc, scalar * tcl, scalar cr, scalar ibm0, 
                     face vector ibmf0, vector ns, scalar alphas, vector nf, scalar alphaf,
                     vector nfh, scalar alphafh, int last)
{
  scalar flux[];
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
  
#if VPRINT
  double cerror1 = 0;
  foreach(reduction(+:cerror1))
    cerror1 += cr[]*dv();
  fprintf (stderr, "initial sweep volume = %0.15g\n", cerror1);
#endif

  scalar ct[];   // adjusted volume of three phase cells
  scalar ctid[]; // corresponding id field

#if _OPENMP
  foreach_face(x, serial) {
#else
  foreach_face(x, reduction (max:cfl)) {
#endif

#if IBM
    double un = uf.x[]*dt/(Delta*fm.x[]*ibmf_temp.x[] + SEPS), s = sign(un);
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
    //coord tempnf = {-s*nf.x[i], nf.y[i], nf.z[i]};
    coord tempnf = {-s*nfh.x[i], nfh.y[i], nfh.z[i]};
    coord lhs = {-0.5, -0.5, -0.5}, rhs = {s*un - 0.5, 0.5, 0.5};

    ct[] = c[];
    ctid[] = 0;

    if (un == 0)
        cf = 0;
    else if (ibm0[i] >= 1.) {
    //else if (ibm0[i]) {
        cf = (cr[i] <= 0. || cr[i] >= 1.)? cr[i] : rectangle_fraction (tempnf, alphafh[i], lhs, rhs);
    }
    #if 1
    else if (ibm0[i] > 0. && ibm0[i] < 1.) {
        coord tempns = {-s*ns.x[i], ns.y[i], ns.z[i]};

        //double advVolume = cfl*ibmf_temp.x[];
        double advVolume = fabs(un)*ibmf_temp.x[];
        //double advVolume = fabs(uf.x[]*dt*pow(Delta, dimension-1)); // note uf.x is already scaled by ibmf & fm
        #if VPRINT && 0
        if (cr[i] && advVolume > cm[i]*ibm0[i]*pow(Delta, dimension))
            fprintf(stderr, "WARNING: (%g, %g, %g) advected area > real area! A_adv = %g A_real = %g\n",
                x, y, z, advVolume, ibm0[i]*pow(Delta, dimension));
        #endif

        //advVolume /= cm[i]*pow(Delta, dimension); // normalize with cell volume

        if (cr[i] <= 0.)
            cf = 0.;
        else if (cr[i] >= ibm0[i]-1e-10) { // interfacial cell is full

            #if VPRINT
            fprintf(stderr, "advected volume_b = %g\n", advVolume);
            #endif

            cf = 1;

            #if VPRINT
            fprintf(stderr, "VOF (b): %g (%g, %g, %g) ibm[%d]=%g cf=%0.15g"
                            " c[%d]=%0.15g cr[%d]=%0.15g un=%g, uf=%g\n", 
                             indicator.x, x, y, z, i, ibm[i], cf, i, c[i], i, cr[i], 
                             un, uf.x[]);
           #endif
        
        }
        else if (cr[i] > 0. && cr[i] < ibm0[i]-1e-10) {
            #if VPRINT
            fprintf(stderr, "advected volume_c = %g\n", advVolume);
            #endif
            coord tempnfh = {-s*nfh.x[i], nfh.y[i], nfh.z[i]};
            
            if (ch[i] >=1 && !tempnfh.x && !tempnfh.y && !tempnfh.z) {
                cf = 1;
            }
            else if (ch[i] <= 1e-10 && !tempnfh.x && !tempnfh.y && !tempnfh.z) {
                cf = 0;
            }
            else {
                double alphacr = immersed_alpha (ch[i], ibm[i], tempnfh, alphafh[i], tempns, alphas[i], cr[i]);
                double newc = plane_volume (tempnfh, alphacr);
                cf = immersed_fraction (newc, tempnfh, alphacr, tempns, alphas[i], lhs, rhs, advVolume, 0);
                
                ct[] = newc;
                ctid[] = i;
            }

            #if VPRINT
            fprintf(stderr, "VOF (c): %g (%g, %g, %g) ibm[%d]=%g cf=%0.15g"
                            " c[%d]=%0.15g cr[%d]=%0.15g un=%0.15g, uf=%g ch=%g\n", 
                             indicator.x, x, y, z, i, ibm[i], cf, i, c[i], i, cr[i], 
                             un, uf.x[], ch[i]);
            #endif // VPRINT
       }
       else
           cf = 0;
    }
    #endif
    else {
        cf = 0;
    }

    flux[] = cf*uf.x[];
    fluxr.x[] = flux[];

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

  boundary({ctid, ct, flux});

  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "src/vof.h:%d: warning: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     __LINE__, cfl - 0.5), fflush (ferr);

  double crsum = 0, crsumclamp = 0;
  foreach(reduction (+:crsum) reduction (+:crsumclamp))
    if (ibm0[] > 0) {

#if VPRINT
      if ((on_interface(ibm) && c[]) || cr[] > ibm0[] || cr[] < 0)
        fprintf (stderr, "F* W/O FLUX: %g (%g, %g, %g) c[]=%0.15g cr[]=%0.15g"
                         " ibmf[]=%g ibmf[1]=%g ibm[]=%g uf[]=%g uf[1]=%g dx=%g\n",
                         indicator.x, x, y, z, c[], cr[], ibmf0.x[], ibmf0.x[1], ibm0[],
                         uf.x[], uf.x[1], Delta);
#endif
      #if 1
      if (ibm0[] > 0. && ibm0[] < 1. && cr[] > 0. && cr[] < ibm0[]-1e-10) // three-phase cell
        for (int _i = -1; _i <= 1; _i += 1) 
            if (ctid[_i] == -_i && ibm0[_i] > 0. && ibm0[_i] < 1.)
                c[] = ct[_i];
      #endif

#if AXI
      double val = cm[];
#else
      double val = 1;
#endif

      c[]  += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[] - divs.x[]))/(val*Delta);
      cr[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[] - divs.x[]))/(val*Delta);

      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
        t[] += ibmCells[]*dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(val*Delta);

#if VPRINT
      if ((on_interface(ibm) && c[]) || cr[] > ibm0[] || cr[] < 0)
        fprintf (stderr, "F* W/FLUX: %g (%g, %g, %g) c[]=%0.15g cr[]=%0.15g cc[]=%g"
                         " flux[]=%g flux[1]=%g div=%g\n",
                         indicator.x, x, y, z, c[], cr[], cc[], flux[], flux[1], divg1[]);
#endif
      crsum += cr[]*pow(Delta, dimension)*val;
      crsumclamp += clamp(cr[], 0, ibm0[])*pow(Delta, dimension)*val;

    }
#if 0
  if (!approx_equal_double (crsum, crsumclamp, 1e-14))
    fprintf (stderr, "WARNING %g: crsum != crsumclamp. crsum=%0.15g crsumclamp=%0.15g err=%g\n", 
        indicator.x, crsum, crsumclamp, get_percent_error(crsum, crsumclamp));
#endif
#if VPRINT
  double cerror2 = 0;
  foreach(reduction(+:cerror2))
    cerror2 += cr[]*dv();
  fprintf (stderr, "(a) post sweep volume = %g cr_sum=%0.15g cr_sumclamp=%0.15g\n", cerror2, crsum, crsumclamp);
#endif

  // update c to conserve cr
#if MOVING
  move_solid_x(ibm0, ibmf0);
  reconstruction_ibm (ibm0, ibmf0, ns, alphas);

  double cerror1a2_x = real_volume(c);

  #if VPRINT
  fprintf (stderr, "(a2) post moving solid volume = %0.15g\n", cerror1a2_x);
  #endif

#endif

  // TODO: only call this function when crsum != crsumclamp!
  double verror = redistribute_volume (c, cr, ibm);
  //double verror = 0;

  #if VPRINT
  fprintf (stderr, "(a3) volume error = %0.15g\n", verror);
  double cerror3 = 0;
  foreach(reduction(+:cerror3))
    cerror3 += cr[]*dv();
  fprintf (stderr, "(b) post sweep volume = %0.15g\n", cerror3);
  #endif
  (void) verror;

  foreach() {
    if (cr[] < 1e-11)
        cr[] = c[] = 0;
    if (on_interface(ibm) && cr[] > ibm[] - 1e-10)
        cr[] = ibm[];
  }

  boundary({c, cr});
  //reconstruction (c, nf, alphaf);

#if VPRINT
  fprintf(stderr, "%g SETTING CONTACT ANGLE TENSION\n", indicator.x);
#endif

  trash({ch});
  foreach() {
    if (ibm[] > 0 && ibm[] < 1 && cr[] >= ibm[]-1e-6)
        ch[] = 1.;
    else
        ch[] = cr[];
    // c[] = cr[];
  }
  boundary({ch});
  reconstruction (ch, nf, alphaf);
  set_contact_angle_tension(ch, cr, ibm0, nf, alphaf, ns, alphas);

  //reconstruction (ch, nf, alphaf);

  if (!last)
      reconstruction (ch, nfh, alphafh);

  delete (tfluxl); free (tfluxl);

}




/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void clean_fluid_real (scalar f, scalar fr, scalar ibm);

void vof_advection (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {
    vector nf[], nfh[], ns[];
    scalar alphaf[], alphafh[], alphas[];

    trash({ch});

    if (i == 0)
        foreach()
            ch[] = c[];

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
    reconstruction_ibm (ibm, ibmf, ns, alphas);
    reconstruction (ch, nfh, alphafh);

    foreach() {
        if (i == 0) {
            if (on_interface(ibm) && on_interface(c)) {
                cr[] = immersed_fraction (c[], (coord){nf.x[], nf.y[], nf.z[]}, alphaf[],
                                         (coord){ns.x[], ns.y[], ns.z[]}, alphas[],
                                         (coord){-0.5, -0.5, -0.5},
                                         (coord){0.5, 0.5, 0.5}, 0) * ibm[];
            }
        else
            cr[] = c[]*ibm[];
        }
        cc[] = (cr[] > 0.5*ibm[]);

        if (on_interface(ibm)) {
#if MOVING
            coord mp, n;
            //double area = ibm0_geometry (point, &mp, &n);
            double area = ibm_geometry (point, &mp, &n);
            double mpx, mpy, mpz;
            local_to_global (point, mp, &mpx, &mpy, &mpz);

            // NOTE: t - dt because uibm_x functions have t + dt, so this is
            // essentially (t - dt) + dt, which cancels out the dt's.
            foreach_dimension()
                divs.x[] = uibm_x(mpx,mpy,mpz, t - dt) * n.x * area;
#else
            foreach_dimension()
                  divs.x[] = 0;
#endif
        }
        else {
            foreach_dimension()
                divs.x[] = 0;
        }
    }

    foreach_face() {
        ibmf_temp.x[] = ibmf.x[];
        uf.x[] *= ibmf_temp.x[];
    }

    boundary({uf});
   
#if VPRINT
    double cerror0 = 0;
    foreach(reduction(+:cerror0))
      cerror0 += cr[]*dv();
    fprintf(stderr, "initial real volume cr = %g\n", cerror0);
#endif

    void (* sweep[dimension]) (scalar, scalar, scalar, scalar *, scalar, scalar, 
                               face vector, vector, scalar, vector, scalar, vector, scalar, int);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++) {

        #if VPRINT
        char ind = (i + d) % dimension == 0? 'x':
                   (i + d) % dimension == 1? 'y': 'z';
        fprintf(stderr, "\n=== %c SWEEP (i = %d) ===\n", ind, i);
        #endif

        int last = (d == dimension - 1);
        sweep[(i + d) % dimension] (c, ch, cc, tcl, cr, ibm, ibmf, ns, alphas, 
                                  nf, alphaf, nfh, alphafh, last);
    }
    delete (tcl), free (tcl);

    foreach_face()
        uf.x[] /= (ibmf_temp.x[] + SEPS);

    boundary({uf});

    //foreach()
    //    c[] = cr[];

    //clean_fluid_real (c, cr, ibm);
#if 0
   reconstruction (c, nfg, alphafg);
   reconstruction (ibm, nsg, alphasg);
#endif
  }
}

event vof (i++)
  vof_advection (interfaces, i);

