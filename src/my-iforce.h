/**
# Interfacial forces

We assume that the interfacial acceleration can be expressed as
$$
\phi\mathbf{n}\delta_s/\rho
$$
with $\mathbf{n}$ the interface normal, $\delta_s$ the interface Dirac
function, $\rho$ the density and $\phi$ a generic scalar field. Using
a CSF/Peskin-like approximation, this can be expressed as
$$
\phi\nabla f/\rho
$$
with $f$ the volume fraction field describing the interface.

The interfacial force potential $\phi$ is associated to each VOF
tracer. This is done easily by adding the following [field
attributes](/Basilisk C#field-attributes). */

attribute {
  scalar phi;
}

/**
Interfacial forces are a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

/**
The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */
scalar phig[];
event acceleration (i++)
{
 

  /**
  We check for all VOF interfaces for which $\phi$ is allocated. The
  corresponding volume fraction fields will be stored in *list*. */

  scalar * list = NULL;
  for (scalar f in interfaces)
    if (f.phi.i) {
      list = list_add (list, f);

      /**
      To avoid undeterminations due to round-off errors, we remove
      values of the volume fraction larger than one or smaller than
      zero. */

      foreach() {
#if CA
	    f[] = clamp (f[], 0., cs[]);
#else
	    f[] = clamp (f[], 0., 1.);
#endif
      }
    }

  /**
  On trees we need to make sure that the volume fraction gradient
  is computed exactly like the pressure gradient. This is necessary to
  ensure well-balancing of the pressure gradient and interfacial force
  term. To do so, we apply the same prolongation to the volume
  fraction field as applied to the pressure field. */
  
#if TREE
  for (scalar f in list) {
    f.prolongation = p.prolongation;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif

  /**
  Finally, for each interface for which $\phi$ is allocated, we
  compute the interfacial force acceleration
  $$
  \phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f
  $$ 
  */

#if 0
  foreach() {
    scalar phi = f.phi;
    if (on_interface(cs) && !f[] && ch[] && !extra[]) {
       // phi[] = 0;
       bool skip = true;
       foreach_neighbor(1) {
         if (extra[])
            skip = false;
       }
       if (skip)
         phi[] = 0;

         if (extra[-1] || extra[1] || extra[0,-1] || extra[0,1]) {
            foreach_dimension() {
                for (int i = -1; i < 2; i += 2)
                    if (extra[i]) {
                        phi[] = phi[i];
                    }
            }
        }
    }
  } 
#endif
    foreach() {
        scalar phi = f.phi;
        phig[] = phi[];
    }

  face vector ia = a;
  foreach_face()
    for (scalar f in list)
      //if (ch[] != ch[-1] && fm.x[] > 0.) {
      if (f[] != f[-1] && fm.x[] > 0.) {

	/**
	We need to compute the potential *phif* on the face, using its
	values at the center of the cell. If both potentials are
	defined, we take the average, otherwise we take a single
	value. If all fails we set the potential to zero: this should
	happen only because of very pathological cases e.g. weird
	boundary conditions for the volume fraction. */

	scalar phi = f.phi;
	double phif =
	  (phi[] < nodata && phi[-1] < nodata) ?
	  (phi[] + phi[-1])/2. :
	  phi[] < nodata ? phi[] :
	  phi[-1] < nodata ? phi[-1] :
	  0.;

    //double val1 = ch[], val2 = ch[-1];
    double val1 = f[]/(cs[] + SEPS), val2 = f[-1]/(cs[-1] + SEPS);
    //double val1 = f[], val2 = f[-1];
    #if 0
    if ((!f[]) && cs[] < 1) {
      bool skip = true;
      double phia = 0;
      int count = 0;
      foreach_neighbor() {
        if (extra[]) {
          skip = false;
        }
        if (f[] && phi[] && extra[]) {
          phia += phi[];
          count++;
        }
      }
      if (skip)
        continue;
     // phif = phia/((double)count);
    }
    #endif

#if IBM || EMBED
    #if CA
      //double val1 = f[]? ch[]: 0;
      //double val2 = f[-1]? ch[-1]: 0;
    #if EMBED
      ia.x[] += alpha.x[]/(fm.x[] + SEPS)*phif*(val1 - val2)/Delta;
    #else
      ia.x[] += alpha.x[]/(fs.x[] + SEPS)*phif*(val1 - val2)/Delta;
    #endif
    #else
      ia.x[] += alpha.x[]/(fs.x[] + SEPS)*phif*(f[] - f[-1])/Delta;
    #endif // CA
#else
      ia.x[] += alpha.x[]/(fm.x[] + SEPS)*phif*(f[] - f[-1])/Delta;
#endif
      }

#if 0
    foreach_face() {
      if (ch[] != ch[-1] && fm.x[] > 0. && fs.x[] < 1 && fs.x[] > 0) {
        if (!f[]) {
            double asum = 0;
            int count = 0;
            foreach_neighbor(1) {
                if (extra[] && fs.x[]) {
                    asum += ia.x[];
                    count++;
                }
            }
            if (count)
                ia.x[] = asum/((double)count);
        }
      }
    }
    #endif

  /**
  On trees, we need to restore the prolongation values for the
  volume fraction field. */
  
#if TREE
  for (scalar f in list) {
    f.prolongation = fraction_refine;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif

  /**
  Finally we free the potential fields and the list of volume
  fractions. */

  for (scalar f in list) {
    scalar phi = f.phi;
    delete ({phi});
    f.phi.i = 0;
  }
  free (list);
}

/**
## References

See Section 3, pages 8-9 of:

~~~bib
@hal{popinet2018, hal-01528255}
~~~
*/
