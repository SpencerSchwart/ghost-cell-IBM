#include "fractions.h"

#define VPRINT 0

// TODO: use ibmf to calculate ibm normal during reconstruction!!!
//      - reconstruction uses MYC, while ibm_geomtry uses ibmf to find n

attribute {
  scalar * tracers, c;
  bool inverse;
}


extern scalar * interfaces;
extern face vector uf;
extern double dt;

double cerror0, cerror1_x, cerror1_y, cerror2, 
       cerror1a_x, cerror1a_y, cerror1b_x, cerror1b_y,
       cerror1c_x, cerror1c_y;

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

void set_contact_angle(scalar f, scalar fr, scalar ibm);

coord indicator = {0,1};

scalar cr[]; // holds the real fluid
scalar cr0[];
scalar c0[];

vector nfg[], nsg[], divs[];
scalar alphafg[], alphasg[];

face vector ibmf_temp[];

scalar cr_x[], cr_y[], cr0_x[], cr0_y[];
vector fluxr[];

void move_solid_x(scalar ibm, face vector ibmf);
void move_solid_y(scalar ibm, face vector ibmf);

foreach_dimension()
static void sweep_x (scalar c, scalar cc, scalar * tcl, scalar cr, scalar ccr, scalar ibm0, face vector ibmf0)
{
  vector nf[], ns[];
  scalar alphaf[], alphas[], flux[];
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
 
  cerror1_x = real_volume (c);

  foreach() {
      cr0_x[] = cr[];
  }

#if VPRINT
  fprintf (stderr, "initial sweep volume = %g\n", cerror1_x);
#endif

  foreach_face(x, reduction (max:cfl)) {

#if IBM
    double un = uf.x[]*dt/(Delta*ibmf_temp.x[] + SEPS), s = sign(un);
    //double un = uf.x[]*dt/(Delta + SEPS), s = sign(un);
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

    if (un == 0)
        cf = 0;
    else if (ibm0[i] >= 1.) {
        cf = (c[i] <= 0. || c[i] >= 1.)? c[i] : rectangle_fraction (tempnf, alphaf[i], lhs, rhs);
    }
    else if (ibm0[i] > 0. && ibm0[i] < 1.) {
        coord tempns = {-s*ns.x[i], ns.y[i], ns.z[i]};
        if (cr[i] <= 0.)
            cf = 0.;
        else if (cr[i] >= ibm0[i]-1e-10) { // interfacial cell is full
            //double alphac = plane_alpha (cr[i], (coord){ns.x[i], ns.y[i]});
            double alphac = alphas[i];
            //cf = (cr[i] <= ibm. || cr[i] >= 1.)? c[i] : rectangle_fraction (tempns, alphac, lhs, rhs);
            //cf = rectangle_fraction (tempns, alphac, lhs, rhs);
            
            coord lhst = {lhs.x,0.5}, rhsb = {rhs.x,-0.5}; 
            coord rect[4] = {lhs,rhsb,rhs,lhst};
            double areaTotal = polygon_area (4, rect); // total area being considered for advection (uf*dt*h)
            double areaLiquid = rectangle_fraction (tempns, alphac, lhs, rhs);

            // TODO: temp fix when flux area > fluid area (small cell problem)
            if (alphac - tempns.x*rhs.x - tempns.y*rhs.y < 0 && alphac - tempns.x*rhsb.x - tempns.y*rhsb.y < 0)
                cf = immersed_fraction (c[i], tempnf, alphaf[i], tempns, alphas[i], lhs, rhs, 0);
            else
                cf = 1;


            #if VPRINT
            fprintf(stderr, "VOF (b): %g (%g, %g) ibm[%d]=%g cf=%0.15g"
                            " c[%d]=%0.15g cr[%d]=%0.15g un=%g, uf=%g\n", 
                             indicator.x, x, y, i, ibm[i], cf, i, c[i], i, cr[i], 
                             un, uf.x[]);
           #endif
        }
        else if (cr[i] > 0. && cr[i] < ibm0[i]) {
            if (ibm0[i] > 0 && ibm0[i] < 1 && cr[i] > 0 && cr[i] < 1) { // three-phase cell
                cf = immersed_fraction (c[i], tempnf, alphaf[i], tempns, alphas[i], lhs, rhs, 0);
            }
            else if (cr[] < 1e-10)
                cf = 0;
            else
                cf = 0;
            #if VPRINT
            fprintf(stderr, "VOF (c): %g (%g, %g) ibm[%d]=%g cf=%0.15g"
                            " c[%d]=%0.15g cr[%d]=%0.15g un=%g, uf=%g\n", 
                             indicator.x, x, y, i, ibm[i], cf, i, c[i], i, cr[i], 
                             un, uf.x[]);
           #endif
       }
       else
           cf = 0;
    }
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
      if (on_interface(ibm) || (ibm0[-1] < 1 && ibm0[-1] > 0) || cr[] > ibm0[] || cr[] < 0)
        fprintf (stderr, "F* W/O FLUX: %g (%g, %g) c[]=%0.15g cr[]=%0.15g"
                         " ibmf[]=%g ibmf[1]=%g ibm[]=%g uf[]=%g uf[1]=%g dx=%g\n",
                         indicator.x, x, y, c[], cr[], ibmf0.x[], ibmf0.x[1], ibm0[],
                         uf.x[], uf.x[1], Delta);
#endif
      c[]  += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[] - divs.x[]))/(Delta);
      cr[] += dt*(flux[] - flux[1] + ccr[]*(uf.x[1] - uf.x[] - divs.x[]))/(Delta);

      cr_x[] = cr[];

#if VPRINT
      if (on_interface(ibm) || (ibm0[-1] < 1 && ibm0[-1] > 0) || cr[] > ibm0[] || cr[] < 0)
        fprintf (stderr, "F* W/FLUX: %g (%g, %g) c[]=%0.15g cr[]=%0.15g cc[]=%g"
                         " flux[]=%g flux[1]=%g div=%g\n",
                         indicator.x, x, y, c[], cr[], cc[], flux[], flux[1], divg1[]);
#endif
      crsum += cr[]*sq(Delta);

      //cr[] = clamp(cr[], 0., ibm0[]);
      crsumclamp += clamp(cr[], 0, ibm0[])*sq(Delta);
    }

  cerror1a_x = real_volume(c);

#if VPRINT
  fprintf (stderr, "(a) post sweep volume = %g cr_sum=%g cr_sumclamp=%g\n", cerror1a_x, crsum, crsumclamp);
#endif

  // update c to conserve cr
  #if MOVING
  move_solid_x(ibm0, ibmf0);

  double cerror1a2_x = real_volume(c);

  #if VPRINT
  fprintf (stderr, "(a2) post moving solid volume = %g\n", cerror1a2_x);
  #endif

  //double verror = redistribute_volume (c, cr, ibm);

  #if VPRINT
  //fprintf (stderr, "(a3) volume error = %0.15g\n", verror);
  #endif

  #endif
  set_contact_angle(c, cr, ibm0);
  reconstruction (c, nf, alphaf);
  reconstruction (ibm0, ns, alphas);

  immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);

  cerror1b_x = real_volume(c);

#if VPRINT
  fprintf (stderr, "(b) post sweep volume = %g\n", cerror1b_x);
#endif

#if 1
  // TODO: have a check here to see if a second reconstruction is necessary
  // TODO: also see if 3rd reconstruction helps at all
  reconstruction (c, nf, alphaf);
  immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);

  cerror1c_x = real_volume(c);

#if VPRINT
  fprintf (stderr, "(c) post sweep volume = %g\n", cerror1c_x);
#endif

#endif

#if 0

  reconstruction (c, nf, alphaf);
  immersed_reconstruction (c, cr, nf, alphaf, ns, alphas);

  double cerror1d_x = real_volume(c);

#if VPRINT
  fprintf (stderr, "(d) post sweep volume = %g\n", cerror1d_x);
#endif

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
                                         (coord){0.5, 0.5, 0.5}, 0) * ibm[];
     }
     else
       creal[] = c[]*ibm[];
     
     cr0[] = creal[];   
     c0[] = c[];
     cc[] = (creal[] > 0.5*ibm[]);
     ccr[] = (creal[] > 0.5*ibm[]);

     if (on_interface(ibm)) {
         coord mp, n;
          //double area = ibm0_geometry (point, &mp, &n);
          double area = ibm_geometry (point, &mp, &n);
          double mpx, mpy, mpz;
          local_to_global (point, mp, &mpx, &mpy, &mpz);

          // NOTE: t - dt because uibm_x functions have t + dt, so this is
          // essentially (t - dt) + dt, which cancels out the dt's.
          divs.x[] = uibm_x(mpx,mpy,mpz, t - dt) * n.x * area;
          divs.y[] = uibm_y(mpx,mpy,mpz, t - dt) * n.y * area;
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
    
    cerror0 = real_volume (c);

    #if VPRINT
    fprintf(stderr, "initial real volume cr = %g\n", cerror0);
    #endif

    void (* sweep[dimension]) (scalar, scalar, scalar *, scalar, scalar, scalar, face vector);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++) {
      char ind = (i + d) % dimension == 0? 'x': 'y';
     // char ind = d == 0? 'x': 'y';

      #if VPRINT
      fprintf(stderr, "\n=== %c SWEEP (i = %d) ===\n", ind, i);
      #endif
      sweep[(i + d) % dimension] (c, cc, tcl, creal, ccr, ibm, ibmf);
      //sweep[d == 1] (c, cc, tcl, creal, ccr, ibm, ibmf);
    }
    delete (tcl), free (tcl);

    foreach_face() {
        uf.x[] /= (ibmf_temp.x[] + SEPS);
    }

   foreach() {
     cr[] = creal[];
   }

   cerror2 = real_volume (c);

    #if VPRINT
    fprintf(stderr, "final 1 real volume cr = %g\n", cerror2);
    #endif

   // TODO: WHY DOES THIS HAVE AN EFFECT ON THE ERROR!?!?!?
   // -update: apparently it doesn't anymore
   #if 0
   reconstruction (c, nf, alphaf);
   reconstruction (ibm, ns, alphas);
   immersed_reconstruction (c, creal, nf, alphaf, ns, alphas);
   #endif

   double cerror2a = real_volume (c);

    #if VPRINT
    fprintf(stderr, "final 2 real volume cr = %g\n", cerror2a);
    #endif

   reconstruction (c, nfg, alphafg);
   reconstruction (ibm, nsg, alphasg);
  }
}

event vof (i++)
  vof_advection (interfaces, i);

