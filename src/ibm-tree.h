/**
### Most of this is taken or inspired from embed-tree.h ###

# Immersed boundaries on adaptive trees

This file defines the restriction/prolongation functions which are
necessary to implement [immersed boundaries](ibm-gcm.h) on adaptive
meshes.

## Volume fraction field *cs*

For the fluid volume fraction field *cs*, the function below is modelled
closely on the volume fraction refinement function
[fraction_refine()](fractions.h#fraction_refine). */

static void ibm_fraction_refine (Point point, scalar cs)
{
  double cc = cs[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      cs[] = cc;
  }
  else {

    /**
    If the cell contains the ibmded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, cs, fs);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	    nc.x = child.x*n.x;
      cs[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}

/**
## Surface fractions field *fs*

The immersed surface fractions *fs* are reconstructed using this
function. */

foreach_dimension()
static void ibm_face_fraction_refine_x (Point point, scalar s)
{
  vector fs = s.v;
  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cs[] <= 0. || cs[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	    fine(fs.x,1,j,k) = cs[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	     (is_local(cell) || is_local(neighbor(2*i-1))))
	    for (int j = 0; j <= 1; j++)
          for (int k = 0; k <= 1; k++)
            fine(fs.x,2*i,j,k) = fs.x[i];
  }
  else {

    /**
    If the cell contains the ibmded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, cs, fs);
    double alpha = plane_alpha (cs[], n);
      
    /**
    We need to reconstruct the face fractions *fs* for the fine cells.
    
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
      fine(fs.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fs.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fs.x,1,0) = fine(fs.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	    if (!fine(cs,0,j,k) || !fine(cs,1,j,k))
	      fine(fs.x,1,j,k) = 0.;
	    else {
	      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
    	  coord nc;
	      nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	      fine(fs.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
	    }

#endif // dimension == 3
    
    /**
    For the fine face fractions coincident with the faces of the
    coarse cell, we obtain the intersection position from the
    coarse cell face fraction. */

    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1).neighbors && (is_local(cell) || is_local(neighbor(2*i-1)))) {
	    if (!is_refined(neighbor(2*i-1))) {
	      if (fs.x[i] <= 0. || fs.x[i] >= 1.)
	        for (int j = 0; j <= 1; j++)
	          for (int k = 0; k <= 1; k++)
		        fine(fs.x,2*i,j,k) = fs.x[i];
	      else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	        double a = fs.y[0,1] <= 0. || fs.y[2*i-1,1] <= 0. ||
	                   fs.y[] >= 1. || fs.y[2*i-1] >= 1.;
	        if ((2.*a - 1)*(fs.x[i] - 0.5) > 0.) {
	          fine(fs.x,2*i,0) = a;
	          fine(fs.x,2*i,1) = 2.*fs.x[i] - a;
	        }
	        else {
	          fine(fs.x,2*i,0) = 2.*fs.x[i] + a - 1.;
    	      fine(fs.x,2*i,1) = 1. - a;
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
		    fine(fs.x,2*i,j,k) =
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
	    if (fine(fs.x,2*i,j,k) && !fine(cs,i,j,k))
	      fine(fs.x,2*i,j,k) = 0.;
      }
  }

}

/**
## Restriction of cell-centered fields

We now define restriction and prolongation functions for cell-centered
fields. The goal is to define second-order operators which do not use
any values from cells entirely contained within the ibmded boundary
(for which *cs = 0*). 

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
  double sum = 0;
  int count = 0;
  foreach_child() {
    sum += s[];
    count += s[]? 1: 0;
  }
  s[] = sum/(count + SEPS);
}


/**
## Refinement/prolongation of cell-centered fields

This is adopted from the embedded boundaries method. For the ghost-cell immersed
boundary method, we want to do something similar, i.e. ignore certain cells when
interpolating near the immersed boundary. However, the cells in which are valid
for interpolation changes. Simply put, any fluid cell or ghost cell
can be used in interpolation. Likewise, we want to avoid using non-ghost solid
cells which have a trivial solution (0). 

Checking to see if a cell is a "fluid" cell is easy: if cs > GCV. However, 
accurately identifying ghost cells is tricker, since cs can equal 0 like solid
cells we want to avoid. The most robust way would be to follow a similar approach
of the is_ghost_cell() function, which requires checks of surrounding cells. 

To avoid this complex implementation, we opt for a simpler method to identify
ghost cells: a valid ghost cell is a cell with cs <= GCV and with a non-trivial
value for the field s (s != 0) of which is getting interpolated.

An interpolation scheme is based on how many bounding neighbors satisfy this
criteria. In 2D, bilinear, triangular, or diagonal interpolation is used if
four, three, or two nodes are fluid or ghost cells. Correspondingly, in 3D, 
trilinear, tetrahedral, and diagonal interpolation is used for 8, 4, and 2
neighbors. */

#define coarse_valid(s,i,j,k) (coarse(cs,i,j,k) > GCV || (coarse(cs,i,j,k) <= GCV && coarse(s,i,j,k)))

static inline void refine_ibm_linear (Point point, scalar s)
{
  foreach_child() {
    if (cs[] < GCV)
      s[] = 0.;
    else {
      assert (coarse(cs));
#if dimension == 2
      if (coarse_valid(s,0,0,0) && coarse_valid(s,child.x,0,0) && 
          coarse_valid(s,0,child.y,0) && coarse_valid(s,child.x,child.y,0)) 
        // bilinear interpolation
        s[] = (9.*coarse(s) + 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + coarse(s,child.x,child.y))/16.;
      else if (coarse_valid(s,0,0,0) && coarse_valid(s,child.x,0,0) && 
               coarse_valid(s,0,child.y,0))
	    // triangular interpolation
        s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      else if (coarse_valid(s,0,0,0) && coarse_valid(s,child.x,child.y,0)) 
        // diagonal interpolation
        s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
#else // dimension == 3
      if (coarse_valid(s,0,0,0) && coarse_valid(s,child.x,0,0) && 
          coarse_valid(s,0,child.y,0) && coarse_valid(s,child.x,child.y,0) && 
          coarse_valid(s,0,0,child.z) && coarse_valid(s,child.x,0,child.z) &&
          coarse_valid(s,0,child.y,child.z) && coarse_valid(s,child.x,child.y,child.z))
        // bilinear/trilinear interpolation
        s[] = (27.*coarse(s) + 
               9.*(coarse(s,child.x) + coarse(s,0,child.y) +
                   coarse(s,0,0,child.z)) + 
               3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
                   coarse(s,0,child.y,child.z)) + 
               coarse(s,child.x,child.y,child.z))/64.;
      else if (coarse_valid(s,0,0,0) && coarse_valid(s,child.x,0,0) && 
               coarse_valid(s,0,child.y,0) && coarse_valid(s,0,0,child.z))
        // tetrahedral interpolation
        s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
               coarse(s,0,0,child.z))/4.;
      else if (coarse_valid(s,0,0,0) && coarse_valid(s,child.x,child.y,child.z)) 
        // diagonal interpolation
        s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
        // Pathological cases, use 1D gradients.
        s[] = coarse(s);
        foreach_dimension() {
          if (coarse(cs,child.x) > GCV || (coarse(cs,child.x) <= GCV && coarse(s,child.x)))
            s[] += (coarse(s,child.x) - coarse(s))/4.;
          else if (coarse(cs,-child.x) > GCV || (coarse(cs,-child.x) <= GCV && coarse(s,-child.x)))
            s[] -= (coarse(s,- child.x) - coarse(s))/4.;
        }
      }
    }
  }
}


foreach_dimension()
static void refine_metric_injection_x (Point point, scalar s)
{
    vector v = s.v;
    double val = on_interface(cs)? 1.: cs[];
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

static void fraction_refine_metric (Point point, scalar s)
{
  double cc = cs[];

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

    coord n = facet_normal (point, cs, fs);
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
