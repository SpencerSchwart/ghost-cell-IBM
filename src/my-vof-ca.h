#include "fractions.h"

attribute {
  scalar * tracers, c;
  bool inverse;
}

scalar ch[];  // Volume fraction field for curvature calculation

scalar c0[];

extern scalar * interfaces;
extern face vector uf;
extern double dt;

foreach_dimension()
static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (cs[-1]) cl /= cs[-1];
  if (cs[])   cc /= cs[];
  if (cs[1])  cr /= cs[1];
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
  ch.refine = ch.prolongation = fraction_refine;
  ch.dirty = true;
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

vector divs[];

const coord indicator = {0,1,2};

#if MOVING
face vector fs_temp[];
#endif

void move_solid_x(scalar cs, face vector fs);
void move_solid_y(scalar cs, face vector fs);
void move_solid_z(scalar cs, face vector fs);

foreach_dimension()
static void sweep_x (scalar c, scalar ch, scalar cc, scalar * tcl, scalar cs0, 
                     face vector fs0, vector ns, scalar alphas, vector nfh, scalar alphafh, int last)
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
  
  foreach_face(x, reduction (max:cfl)) {

    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
    if (cs[] >= 1.)
#elif IBM
    if (cs0[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    double cf = 0;
    coord tempnf  = {-s*nfh.x[i], nfh.y[i], nfh.z[i]};
    coord lhs = {-0.5, -0.5, -0.5}, rhs = {s*un - 0.5, 0.5, 0.5};

    if (un == 0)
        cf = 0;
    else if (cs0[i] >= 1.) {
        cf = (c[i] <= 0. || c[i] >= 1.)? c[i] : rectangle_fraction (tempnf, alphafh[i], lhs, rhs);
    }
    else if (cs0[i] > 0. && cs0[i] < 1.) {
        coord tempns = {-s*ns.x[i], ns.y[i], ns.z[i]};

        double advVolume = fabs(un)*fs.x[];

        if (c[i] <= 0.)
            cf = 0.;
        else if (c[i] >= cs0[i]-INT_TOL) // interfacial cell is full
            cf = 1;
        else if (c[i] > 0. && c[i] < cs0[i]-INT_TOL) { // three-phase
            double alpha = alphafh[i];
            if (!tempnf.x && !tempnf.y && !tempnf.z) {
                coord off = {0,0,0};
                off.x = i;
                tempnf = youngs_normal_off (point, c, off);
                tempnf = (coord){-s*tempnf.x, tempnf.y, tempnf.z};
                alpha = 0; // inital guess
            }
            
            normalize2(&tempns);
            coord nc = normal_contact (tempns, tempnf, contact_angle[i]);
            normalize_sum(&tempns);
            normalize_sum(&nc);

            tempnf = nc;
            
            double alphacr = immersed_alpha (ch[i], cs[i], tempnf, alpha, tempns, alphas[i], c[i]);
            double newc = plane_volume (tempnf, alphacr);
            cf = immersed_fraction (newc, tempnf, alphacr, tempns, alphas[i], lhs, rhs, advVolume, 0);
       }
       else
           cf = 0;
    }
    else 
        cf = 0;

#if IBM
    flux[] = cf*uf.x[]*fs.x[];
#else
    flux[] = cf*uf.x[];
#endif

    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = cs[i]? c[i]/cs[i]: 0;
      if (t.inverse)
    	cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	    double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
#if IBM
    	tflux[] = ff*cf1*uf.x[]*fs.x[];
#else
     	tflux[] = ff*cf1*uf.x[];       
#endif
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

  double crsum = 0, crsum_clamp = 0;
  foreach(reduction (+:crsum) reduction (+:crsum_clamp)) {

#if AXI
    if (cs0[] > 0 && cm[]) {
      double val = cm[];
      #if EMBED
      val /= cs0[];
      #endif
#else
    if (cs0[] > 0) {
      double val = 1;
#endif

#if RAIN
      if (empty_interfacial(point, c))
      {
        /**
        interfacial cell must border a full cell to avoid creating whisps*/
        bool real_interface = false;
        foreach_neighbor(1) // only 3x3
            if (c[] >= cs[]) {
                real_interface = true; 
                break;
            }

        if (real_interface) {   
            coord n = interface_normal (point, c), p;
            double alpha = plane_alpha (c[], n);
            plane_area_center (n, alpha, &p);
            coord pc = {x, y, z};
            foreach_dimension() // put the midpoint in global coordinates
                pc.x += p.x*Delta;

            double frain = -c.rain.u.x * c.rain.rvf * rain_probability(point, pc, c);

            c[]  += dt*frain/(val*Delta);
            ch[] += dt*frain/(val*Delta);
        }
      }
#endif

#if IBM
      c[]  += dt*(flux[] - flux[1] + cc[]*(fs.x[1]*uf.x[1] - fs.x[]*uf.x[] - divs.x[]))/(val*Delta);
      ch[] += dt*(flux[] - flux[1] + cc[]*(fs.x[1]*uf.x[1] - fs.x[]*uf.x[] - divs.x[]))/(val*Delta);

      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
        t[] += dt*(tflux[] - tflux[1] + tc[]*(fs.x[1]*uf.x[1] - fs.x[]*uf.x[]))/(val*Delta);
#elif EMBED
      c[]  += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[] - divs.x[]))/(val*Delta);
      ch[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[] - divs.x[]))/(val*Delta);   

      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
        t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(val*Delta);
#endif

      crsum += c[]*pow(Delta, dimension)*val;
      crsum_clamp += clamp(c[], 0., 1.)*pow(Delta, dimension)*val;
    }
  }

#if MOVING
  move_solid_x(cs0, fs0);
  reconstruction_cs (cs0, fs0, ns, alphas);
#endif

  if (crsum != crsum_clamp) {
    fprintf(stderr, "crsum != crsum_clamp: %0.14g != %0.14g\nRedistributing volume...\n",
            crsum, crsum_clamp);
    //redistribute_volume(c, cs);
    foreach() {
        if (c[] > cs[] || c[] < 0)
            c[] = clamp(c[], 0, cs[]);
    }
  }

  if (!last) {
      reconstruction (ch, nfh, alphafh);
  }
  else {
    
    if (c.wetting.dynamic) // if using dynamic CA model, update CA in three-phase cells
        update_contact_angle (c, c0, cs0, ns, alphas, u, contact_angle);

    foreach() {
      if (c[] < VTOL)
        c[] = 0;
      if (cs[] > 0 && cs[] < 1 && c[] >= cs[]-INT_TOL)
          ch[] = 1;
      else if (cs[] > 0 && c[] < 1 && c[])
          ch[] = c[]/cs[];
      else if (cs[] > 0)
          ch[] = c[];
      else
          ch[] = 0;
    }

    scalar alphaf[];
    vector nf[];

    reconstruction(ch, nf, alphaf); // should be based on c or ch? TODO: is it even necessary?
    set_contact_angle(ch, c, cs0, nf, alphaf, ns, alphas); // find ch
  }

  delete (tfluxl); free (tfluxl);
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */

void vof_advection (scalar * interfaces, int i)
{
  for (scalar c in interfaces) {
    vector nf[], nfh[], ns[];
    scalar alphaf[], alphafh[], alphas[];

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
	      tc[] = c[] < 0.5*cs[] && cs[]? t[]/(1. - c[]/cs[]) : 0.;
	    else
	      tc[] = c[] > 0.5*cs[] && cs[] ? t[]/(c[]/cs[]) : 0.;
      }
#endif // !NO_1D_COMPRESSION
    }

    reconstruction_cs (cs, fs, ns, alphas);
    reconstruction (ch, nfh, alphafh);
    if (i == 0) { // is this necessary still?
        reconstruction (c, nf, alphaf);
        real_fluid (ch, c, nf, alphaf, ns, alphas);
    }

    foreach() {
        c0[] = c[];
        cc[] = (c[] > 0.5*cs[]);

        foreach_dimension()
            divs.x[] = 0;
    }

#if MOVING
    foreach_face() {
        fs_temp.x[] = fs.x[];
        uf.x[] *= fs_temp.x[];
    }
    boundary({uf.x, uf.y, uf.z});
#endif

    void (* sweep[dimension]) (scalar, scalar, scalar, scalar *, scalar, 
                               face vector, vector, scalar, vector, scalar, int);
    int d = 0;
    foreach_dimension()
      sweep[d++] = sweep_x;
    for (d = 0; d < dimension; d++) {
        int last = (d == dimension - 1);
        sweep[(i + d) % dimension] (c, ch, cc, tcl, cs, fs, ns, alphas, nfh, alphafh, last);
    }
    delete (tcl), free (tcl);

#if MOVING
    foreach_face()
        uf.x[] /= (fs_temp.x[] + SEPS);
    boundary({uf.x, uf.y, uf.z, c});
#endif
  }
}

event vof (i++)
  vof_advection (interfaces, i);

