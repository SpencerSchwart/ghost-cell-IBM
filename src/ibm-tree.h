/**

###### Most of this is taken from embed-tree.h ######

# Embedded boundaries on adaptive trees

This file defines the restriction/prolongation functions which are
necessary to implement [ibmded boundaries](ibm.h) on adaptive
meshes.

## Volume fraction field *ibm*

For the ibmded fraction field *ibm*, the function below is modelled
closely on the volume fraction refinement function
[fraction_refine()](fractions.h#fraction_refine). */

static void ibm_fraction_refine (Point point, scalar ibm)
{
  double cc = ibm[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      ibm[] = cc;
  }
  else {

    /**
    If the cell contains the ibmded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, ibm, ibmf);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	    nc.x = child.x*n.x;
      ibm[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}

/**
## Surface fractions field *ibmf*

The ibmded surface fractions *ibmf* are reconstructed using this
function. */

foreach_dimension()
static void ibm_face_fraction_refine_x (Point point, scalar s)
{
  vector ibmf = s.v;

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (ibm[] <= 0. || ibm[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	    fine(ibmf.x,1,j,k) = ibm[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	     (is_local(cell) || is_local(neighbor(2*i-1))))
	    for (int j = 0; j <= 1; j++)
          for (int k = 0; k <= 1; k++)
            fine(ibmf.x,2*i,j,k) = ibmf.x[i];
  }
  else {

    /**
    If the cell contains the ibmded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, ibm, ibmf);
    double alpha = plane_alpha (ibm[], n);
      
    /**
    We need to reconstruct the face fractions *ibmf* for the fine cells.
    
    For the fine face fractions contained within the coarse cell,
    we compute the intersections directly using the VOF
    reconstruction. */

#if dimension == 2

    /**
    In 2D, we obtain the face fractions by taking into
    account the orientation of the normal. */

    if (2.*fabs(alpha) < fabs(n.y)) {
      double yc = alpha/n.y;
      int i = yc > 0.;
      fine(ibmf.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(ibmf.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(ibmf.x,1,0) = fine(ibmf.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	    if (!fine(ibm,0,j,k) || !fine(ibm,1,j,k))
	      fine(ibmf.x,1,j,k) = 0.;
	    else {
	      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
    	  coord nc;
	      nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	      fine(ibmf.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
	    }

#endif // dimension == 3
    
    /**
    For the fine face fractions coincident with the faces of the
    coarse cell, we obtain the intersection position from the
    coarse cell face fraction. */

    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1).neighbors &&
	      (is_local(cell) || is_local(neighbor(2*i-1)))) {
	    if (!is_refined(neighbor(2*i-1))) {
	      if (ibmf.x[i] <= 0. || ibmf.x[i] >= 1.)
	        for (int j = 0; j <= 1; j++)
	          for (int k = 0; k <= 1; k++)
		        fine(ibmf.x,2*i,j,k) = ibmf.x[i];
	      else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	        double a = ibmf.y[0,1] <= 0. || ibmf.y[2*i-1,1] <= 0. ||
	                   ibmf.y[] >= 1. || ibmf.y[2*i-1] >= 1.;
	        if ((2.*a - 1)*(ibmf.x[i] - 0.5) > 0.) {
	          fine(ibmf.x,2*i,0) = a;
	          fine(ibmf.x,2*i,1) = 2.*ibmf.x[i] - a;
	        }
	        else {
	          fine(ibmf.x,2*i,0) = 2.*ibmf.x[i] + a - 1.;
    	      fine(ibmf.x,2*i,1) = 1. - a;
	        }

#else  // dimension == 3

	    /**
	    In 3D we reconstruct the face fraction from the projection
	    of the cell interface reconstruction, as above. */
	  
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++) {
		    static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
		    coord nc;
		    nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
		    fine(ibmf.x,2*i,j,k) =
		      rectangle_fraction (nc, alpha - n.x*(2.*i - 1.)/2., a, b);
	      }

#endif // dimension == 3
	  }
	}

	/**
	The face fractions of empty children cells must be zero. */
	
	for (int j = 0; j <= 1; j++)
	#if dimension > 2
	  for (int k = 0; k <= 1; k++)
	#endif
	    if (fine(ibmf.x,2*i,j,k) && !fine(ibm,i,j,k))
	      fine(ibmf.x,2*i,j,k) = 0.;
      }
  }
}

/**
## Restriction of cell-centered fields

We now define restriction and prolongation functions for cell-centered
fields. The goal is to define second-order operators which do not use
any values from cells entirely contained within the ibmded boundary
(for which *ibm = 0*). 

When restricting it is unfortunately not always possible to obtain a
second-order interpolation. This happens when the parent cell does not
contain enough child cells not entirely contained within the ibmded
boundary. In these cases, some external information (i.e. a boundary
gradient condition) is required to be able to maintain second-order
accuracy. This information can be passed by defining the
*ibm_gradient()* function of the field being restricted. */

attribute {
  void (* ibm_gradient) (Point, scalar, coord *);
}

static inline void restriction_ibm_linear (Point point, scalar s)
{  
  // 0 children
  if (!ibm[]) {
    s[] = 0.;
    return;
  }

  /**
  We first try to interpolate "diagonally". If enough child cells are
  defined (i.e. have non-zero ibmded fractions), we return the
  corresponding value. */

  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)
#if dimension > 2
    for (int j = 0; j <= 1; j++)
#endif
      if (fine(ibm,0,i,j) && fine(ibm,1,!i,!j))
	val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    s[] = val/nv;
    return;
  }

  /**
  Otherwise, we use the average of the child cells which are defined
  (there is at least one). */
  
  coord p = {0.,0.,0.};
  foreach_child()
    // if (ibm[] || is_ghost_cell (point, ibm)) // does including ghost cells help?
    if (ibm[])
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  assert (nv > 0.);
  s[] = val/nv;

  /**
  If the gradient is defined and if the variable is not using
  homogeneous boundary conditions, we improve the interpolation using
  this information. */
  
  if (s.ibm_gradient && s.boundary[0] != s.boundary_homogeneous[0]) {
    coord o = {x,y,z}, g;
    s.ibm_gradient (point, s, &g);
    foreach_dimension()
      s[] += (o.x - p.x/nv)*g.x;
  }
}

/**
## Refinement/prolongation of cell-centered fields

For refinement, we use either bilinear interpolation, if the required
four coarse cell values are defined or trilinear interpolation if only
three coarse cell values are defined. If less than three coarse cell
values are defined ("pathological cases" below), we try to estimate
gradients in each direction and add the corresponding correction. */

static inline void refine_ibm_linear (Point point, scalar s)
{
  foreach_child() {
    if (!ibm[])
      s[] = 0.;
    else {
      assert (coarse(ibm));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(ibmf.x,i) && coarse(ibmf.y,0,j) &&
	  (coarse(ibm) == 1. || coarse(ibm,child.x) == 1. ||
	   coarse(ibm,0,child.y) == 1. || coarse(ibm,child.x,child.y) == 1.)) {
	assert (coarse(ibm,child.x) && coarse(ibm,0,child.y));
	if (coarse(ibmf.x,i,child.y) && coarse(ibmf.y,child.x,j)) {
	  // bilinear interpolation
	  assert (coarse(ibm,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(ibm,child.x,child.y) &&
	       ((coarse(ibmf.x,i) && coarse(ibmf.y,child.x,j)) ||
		(coarse(ibmf.y,0,j) && coarse(ibmf.x,i,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(ibmf.x,i) > 0.25 && coarse(ibmf.y,0,j) > 0.25 &&
	  coarse(ibmf.z,0,0,k) > 0.25 &&
	  (coarse(ibm) == 1. || coarse(ibm,child.x) == 1. ||
	   coarse(ibm,0,child.y) == 1. || coarse(ibm,child.x,child.y) == 1. ||
	   coarse(ibm,0,0,child.z) == 1. || coarse(ibm,child.x,0,child.z) == 1. ||
	   coarse(ibm,0,child.y,child.z) == 1. ||
	   coarse(ibm,child.x,child.y,child.z) == 1.)) {
	assert (coarse(ibm,child.x) && coarse(ibm,0,child.y) &&
		coarse(ibm,0,0,child.z));
	if (coarse(ibmf.x,i,child.y) && coarse(ibmf.y,child.x,j) &&
	    coarse(ibmf.z,child.x,child.y,k) &&
	    coarse(ibmf.z,child.x,0,k) && coarse(ibmf.z,0,child.y,k)) {
	  assert (coarse(ibm,child.x,child.y) && coarse(ibm,child.x,0,child.z) &&
		  coarse(ibm,0,child.y,child.z) &&
		  coarse(ibm,child.x,child.y,child.z));
	  // bilinear interpolation
	  s[] = (27.*coarse(s) + 
		 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		     coarse(s,0,0,child.z)) + 
		 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		     coarse(s,0,child.y,child.z)) + 
		 coarse(s,child.x,child.y,child.z))/64.;
	}
	else
	  // tetrahedral interpolation
	  s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
		 coarse(s,0,0,child.z))/4.;
      }
      else if (coarse(ibm,child.x,child.y,child.z) &&
	       ((coarse(ibmf.z,child.x,child.y,k) &&
		 ((coarse(ibmf.x,i) && coarse(ibmf.y,child.x,j)) ||
		  (coarse(ibmf.y,0,j) && coarse(ibmf.x,i,child.y))))
		||
		(coarse(ibmf.z,0,0,k) &&
		 ((coarse(ibmf.x,i,0,child.z) && coarse(ibmf.y,child.x,j,child.z)) ||
		  (coarse(ibmf.y,0,j,child.z) && coarse(ibmf.x,i,child.y,child.z))))
		||
		(coarse(ibmf.z,child.x,0,k) &&
		 coarse(ibmf.x,i) && coarse(ibmf.y,child.x,j,child.z))
		||
		(coarse(ibmf.z,0,child.y,k) &&
		 coarse(ibmf.y,0,j) && coarse(ibmf.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(ibmf.x,(child.x + 1)/2) && coarse(ibm,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(ibmf.x,(- child.x + 1)/2) && coarse(ibm,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  }
}

/**
## Refinement/prolongation of face-centered velocity

This function is modelled on
[*refine_face_x()*](/src/grid/tree-common.h#refine_face_x) and is
typically used to refine the values of the face-centered velocity
field *uf*. It uses linear interpolation, taking into account the
weighting by the ibmded fractions *ibmf*. */

foreach_dimension()
void refine_ibm_face_x (Point point, scalar s)
{
  vector v = s.v;
  for (int i = 0; i <= 1; i++)
    if (neighbor(2*i - 1).neighbors &&
	(is_local(cell) || is_local(neighbor(2*i - 1)))) {
      double g1 = ibmf.x[i] >= 1. && ibmf.x[i,+1] && ibmf.x[i,-1] ?
	(v.x[i,+1]/ibmf.x[i,+1] - v.x[i,-1]/ibmf.x[i,-1])/8. : 0.;
      double g2 = ibmf.x[i] >= 1. && ibmf.x[i,0,+1] && ibmf.x[i,0,-1] ?
	(v.x[i,0,+1]/ibmf.x[i,0,+1] - v.x[i,0,-1]/ibmf.x[i,0,-1])/8. : 0.;
      for (int j = 0; j <= 1; j++)
	for (int k = 0; k <= 1; k++)
	  fine(v.x,2*i,j,k) = ibmf.x[i] ?
	    fine(ibmf.x,2*i,j,k)*(v.x[i]/ibmf.x[i] +
				(2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
    }
  if (is_local(cell)) {
    double g1 = (ibmf.x[0,+1] + ibmf.x[1,+1]) && (ibmf.x[0,-1] + ibmf.x[1,-1]) ?
      ((v.x[0,+1] + v.x[1,+1])/(ibmf.x[0,+1] + ibmf.x[1,+1]) -
       (v.x[0,-1] + v.x[1,-1])/(ibmf.x[0,-1] + ibmf.x[1,-1]))/8. : 0.;
    double g2 = (ibmf.x[1,0,+1] + ibmf.x[0,0,+1]) && (ibmf.x[1,0,-1] + ibmf.x[0,0,-1]) ?
      ((v.x[0,0,+1] + v.x[1,0,+1])/(ibmf.x[1,0,+1] + ibmf.x[0,0,+1]) -
       (v.x[0,0,-1] + v.x[1,0,-1])/(ibmf.x[1,0,-1] + ibmf.x[0,0,-1]))/8. : 0.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(v.x,1,j,k) = ibmf.x[] + ibmf.x[1] ?
	  fine(ibmf.x,1,j,k)*((v.x[] + v.x[1])/(ibmf.x[] + ibmf.x[1]) +
			    (2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
  }
}

foreach_dimension()
static void refine_metric_injection_x (Point point, scalar s)
{
    vector v = s.v;
    double val = on_interface(ibm)? 1.: ibm[];
    foreach_child()
        v.x[] = val;
}


static inline void face_max_metric (Point point, vector v)
{
  foreach_dimension() {
    #if dimension == 2
      v.x[] = max(fine(v.x,0,0), fine(v.x,0,1));
      v.x[1] = max(fine(v.x,2,0), fine(v.x,2,1));
    #else // dimension == 3
      v.x[] =  max(max(fine(v.x,0,0,0), fine(v.x,0,1,0)),
	               max(fine(v.x,0,0,1), fine(v.x,0,1,1)));
      v.x[1] = max(max(fine(v.x,2,0,0), fine(v.x,2,1,0)),
                   max(fine(v.x,2,0,1), fine(v.x,2,1,1)));
    #endif
  }
}


static inline void restriction_face_metric (Point point, scalar s)
{
  face_max_metric (point, s.v);
}


static inline void restriction_cell_metric (Point point, scalar s)
{
    double sum = 0.;
    foreach_child()
        sum += ibm[];
    s[] = sum/(1 << dimension) > 0.5;
}


static void fraction_refine_metric (Point point, scalar s)
{
  double cc = ibm[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      s[] = cc;
  }
  else {

    /**
    If the cell contains the ibmded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, ibm, ibmf);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	    nc.x = child.x*n.x;
      s[] = rectangle_fraction (nc, alpha, a, b) > 0.5;
    }
  }
}

