/**
# Bell-Collela-Glaz advection scheme

The function below implements the 2nd-order, unsplit, upwind scheme of
[Bell-Collela-Glaz, 1989](references.bib#bell89). Given a centered
scalar field *f*, a face vector field *uf* (possibly weighted by a
face metric), a timestep *dt* and a source term field *src*, it fills
the face vector field *flux* with the components of the advection
fluxes of *f*. */

void tracer_fluxes (scalar f,
		    face vector uf,
		    face vector flux,
		    double dt,
		    (const) scalar src)
{

  /**
  We first compute the cell-centered gradient of *f* in a locally-allocated
  vector field. */
  
  vector g[];

/*
###### GRADIENTS MUST BE CHANGED ######
*/
  gradients ({f}, {g});

  /**
  For each face, the flux is composed of two parts... */

  coord indicator = {1,2};
  foreach_face() {

    /**
    A normal component... (Note that we cheat a bit here, `un` should
    strictly be `dt*(uf.x[i] + uf.x[i+1])/((fm.x[] +
    fm.x[i+1])*Delta)` but this causes trouble with boundary
    conditions (when using narrow '1 ghost cell' stencils)). */
#if IBM
    double metric = fm.x[]? 1: 0;
    double un = metric*dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
#else
    double un = dt*uf.x[]/(fm.x[]*Delta + SEPS), s = sign(un);
#endif

    int i = -(s + 1.)/2.;
    double f2 = f[i] + (src[] + src[-1])*dt/4. + s*(1. - s*un)*g.x[i]*Delta/2.;

    /**
    and tangential components... */

    #if dimension > 1
    if (fm.y[i] && fm.y[i,1]) {
      double vn = (uf.y[i] + uf.y[i,1])/(fm.y[i] + fm.y[i,1]);
      double fyy = vn < 0. ? f[i,1] - f[i] : f[i] - f[i,-1];
      f2 -= dt*vn*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i] && fm.z[i,0,1]) {
      double wn = (uf.z[i] + uf.z[i,0,1])/(fm.z[i] + fm.z[i,0,1]);
      double fzz = wn < 0. ? f[i,0,1] - f[i] : f[i] - f[i,0,-1];
      f2 -= dt*wn*fzz/(2.*Delta);
    }
    #endif
    flux.x[] = f2*uf.x[];
#if 0
    //if (on_interface (ibm))
        fprintf (stderr, "%g %g %g %g %g %g %d f2=%g un=%g g.x[i]=%g u.x[i]=%g uf.x[]=%g fm[i]=%g fm[i,1]=%g flux=%g\n",
                          indicator.x, x, y, ibm[], fm.x[], fm.y[],
                          level, f2, un, g.x[i], f[i], uf.x[], fm.y[i], fm.y[i,1], flux.x[]);
#endif
  }


}

/**
The function below uses the *tracer_fluxes* function to integrate the
advection equation, using an explicit scheme with timestep *dt*, for
each tracer in the list. */

void advection (scalar * tracers, face vector u, double dt,
		scalar * src = NULL)
{

  /**
  If *src* is not provided we set all the source terms to zero. */
  
  scalar * psrc = src;
  if (!src)
    for (scalar s in tracers) {
      const scalar zero[] = 0.;
      src = list_append (src, zero);
    }
  assert (list_len (tracers) == list_len (src));

  scalar f, source;
  for (f,source in tracers,src) {
    face vector flux[];
    #if 0
    foreach_dimension() {
      flux.x.refine = flux.x.prolongation = refine_ibm_face_x;
      flux.x.restriction = restriction_ibm_linear;
    }
    restriction({flux});
    #endif
    tracer_fluxes (f, u, flux, dt, source);
    //boundary({flux});
    //coord indicator = {1,2};
#if !EMBED
    foreach()
      foreach_dimension() {
        // note f[] is u.x[] or u.y[]
#if IBM
#if 0
        if (on_interface(ibm))
            fprintf(stderr, "before %g %g %g %g ibmF.x%g ibmF.y=%g %d u=%g\n", 
                              indicator.x, x, y,ibm[], fm.x[], fm.y[], level, f[]);
#endif                              
        f[] += cm[]*dt*(flux.x[] - flux.x[1])/(Delta*cm[]+SEPS);

        if (fabs(f[]) > LIMIT)
            fprintf(stderr, "WARNING in bcg.h: f[] = %g in (%g, %g) exceeds %g\n", f[], x, y, LIMIT);
#if 0       
        if (on_interface(ibm))
            fprintf(stderr, "after %g %g %g %g ibmF.x%g ibmF.y=%g %d u=%g f[]=%g f[1]=%g fm[1]=%g\n", 
                              indicator.x, x, y,ibm[], fm.x[], fm.y[], level, f[], flux.x[], flux.x[1], fm.x[1]);
#endif
#else
        f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
#endif // IBM
      }
#else // EMBED
    update_tracer (f, u, flux, dt);
#endif // EMBED
  }

  if (!psrc)
    free (src);
}
