#if !AXI
#undef dv
#define dv()  (pow(Delta,dimension))
#endif // !AXI

#include "fractions.h"
#include "ibm-utils.h"
#include "mls.h"

#define BGHOSTS 2
#undef SEPS
#define SEPS 1e-30

#define IBM 1

/**
The GCV variable is the minimum fluid volume before a cell is classified as a
ghost cell. The standard GC-IBM corresponds to GCV = 0.5. In otherwords, if
cs > GCV, then it is a fluid cell. Otherwise, it is a ghost or solid cell.

For very small cells, cs can be << 1. Small cells can cause instability issues
if the momentum equation is solved in them, which is done in the present IBM.
Therefore, GCV is utilized as an emperical cutoff. 1e-4 or 1e-3 seems to work well.

User can also define GCV themselves in the source file */

#ifndef GCV
#define GCV 1e-4 
#endif


/**
# Volume fraction fields

Scalar for storing the solid geomery and corresponding face vector fields. If the
immersed boundary is moving, cs0 is used to store the position of the solid
interface at the previous timestep. 

cs is the fluid volume fraction, i.e.
cs = 1 is a pure fluid cell,
cs = 0 is a pure solid cell,
0 < cs < 1 is an interfacial/cut cell. */

scalar cs[], cs0[];
face vector fs[], fs0[];

/**
The IBM can be used as a symmetric plane to reduce the size of the simulation 
domain. The notation for csp is similiar to cs:

csp = 1 is a fluid cell
csp = 0 is a solid cell. 

The symmetric plane must be defined so that it perfectly coincides with the cell
faces in the domain. When using AMR, the coordinate of the plane must be a
multiple of the largest cell size. */

scalar csp[];

/**
# Metric fields

Cells with gc = 1 are fluid cells, where the governing equations are solved.
Cells with gc = 0 are ghost or solids cells, whose values are either trivial or
determined using interpolation. */

scalar gc[];
face vector gcf[];

double (* metric_ibm_factor) (Point, coord) = NULL; // for axi

typedef struct PointIBM {
    int i, j, k;
} PointIBM;

/**
We define the names that are used to set boundary conditions on the immersed
and symmetric cut plane, corresponding to cs and csp, respectively. */

bid immersed;
bid cut;

attribute {
    vector mp;     // boundary intercepts (TODO: should name this better)

    struct ibm {
        coord (*normal)(coord);
        coord (*bi)(coord);
    } ibm;
}

/**
The bondary condition macros are overridden to allow compatability with the IBM.
The macros make use of ibm_area_center to retrieve the boundary intercept from
ghost cells. */

static inline
double ibm_area_center (Point point, scalar s, double* x1, double* y1, double* z1)
{
    vector mp = s.mp;
    *x1 = mp.x[], *y1 = mp.y[], *z1 = mp.z[];
    return 1;
}

macro2
double dirichlet (double expr, Point point = point,
                  scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = true, expr : 2.*expr - s[];
}

macro2
double dirichlet_homogeneous (double expr, Point point = point,
                              scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = true, 0 : -s[];
}

macro2
double neumann (double expr, Point point = point,
                scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = false, expr : Delta*expr + s[];
}

macro2
double neumann_homogeneous (double expr, Point point = point,
                            scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = false, 0 : s[];
}

macro2
double navier (double expr, Point point = point,
               scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = true, ((bool *)data)[1] = true, expr : 2.*expr - s[];
}

macro2
double navier_homogeneous (double expr, Point point = point,
                           scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = true, ((bool *)data)[1] = true, expr : 2.*expr - s[];
}

#if 0
// Working in progress
macro2
double contact_angle (double expr, Point point = point,
                      scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = true, expr : expr;
}
#endif

/**
This function returns true if the given point has a direct neighbor (N, S, E, W)
that has no liquid volume fraction, i.e. cs == 0, and fills pc and n with the
midpoint and corresponding normal, respectively. 

This is only used if the immersed boundary perfectly coincides with a cell face. */

bool empty_neighbor (Point point, coord * pc, coord * n, scalar cs)
{
    coord pc_temp = {x, y, z}, cc = {x, y, z};
    double cs_temp = cs[];
    double max_d = 1e6;
    int neighbor = 0;

    foreach_direct_neighbor() {
        double distance = distance3D(x - cc.x, y - cc.y, z - cc.z);
        if (cs[] == 0 && cs_temp == 1 && distance < max_d) {
            pc_temp.x = 0.5*(cc.x + x);
            pc_temp.y = 0.5*(cc.y + y);
            pc_temp.z = 0.5*(cc.z + z);

            max_d = distance3D(x - cc.x, y - cc.y, z - cc.z);

            neighbor = 1;
            *pc = pc_temp;
        }
    }
   
    /**
    Calculate the normal. n can only be {1,0}, {-1,0}, {0,1}, or {0,-1}.
    Should probably not compare floating point values like this. */

    coord n_temp = {pc_temp.x != cc.x, 
                    pc_temp.y != cc.y,
                    pc_temp.z != cc.z};
    foreach_dimension() {
        n_temp.x *= sign2(pc_temp.x - cc.x);
    }

    *n = n_temp;

    return neighbor;
}


/**
Checks to see if the given point has at least 1 (direct) fluid neighbor touching
one of it's faces. If so, returns true. */

bool fluid_neighbor (Point point, scalar cs)
{
    foreach_direct_neighbor()
        if (cs[] > GCV && csp[] > GCV)
            return true;
    return false;
}


/**
Using AMR slightly complicates the definition of ghost cells. A naive treatment
will cause coarse cells, behind a layer of finer ghost cells, to be marked as
ghost cells, which are useless and cause problems with mesh constant adaptation.
This is because the foreach_neighbor macro(s) iterate over the cells in the 
stencil with the same level as the original one, using prolongation/restriction 
functions to obtain the values if they are not leaf cells.

Considering the above, we add additional checks in the ghost cell definition,
taken care of by this function.

TODO: probably only need one function, not both fluid_neighbor() and 
      match_level(). */

bool match_level (Point point, scalar cs)
{
#if TREE
    foreach_direct_neighbor() {
        
        /**
        If a neighboring fluid cell is a leaf or coarser than itself, then
        that point is valid and the cell is a ghost cell. */
        if (cs[] > GCV && (is_leaf(cell) || is_prolongation(cell)))
            return true;

        /**
        Otherwise, if the neighboring fluid cell is more refined than itself,
        i.e., has children that are leaves, then we check to see if any of the
        cells are already ghost cells. If less than two are, then we return
        true. */
        else if (cs[] > GCV && is_refined(cell)) {
            
            int numgc = 0;
            for (int i = 0; i <= 1; i++)
                for (int j = 0; j <= 1; j++)
#if dimension == 2
                    if (fine(cs,i,j,0) < GCV)
                        numgc++;
#else
                    for (int k = 0; k <= 1; k++)
                        if (fine(cs,i,j,k) < GCV)
                            numgc++;
#endif                   

            if (numgc < 2)
                return true;
        }
    }
    return false;
#else
    
    /**
    If not using AMR/Trees, then there is no need for this additional check. */
    return true;
#endif
}


/**
is_ghost_cell returns true if the given cell shares a face with a fluid cell,
if cs <= GCV, and if the neighboring fluid cell is valid when using AMR. */

bool is_ghost_cell (Point point, scalar cs)
{
    return cs[] <= GCV && csp[] > GCV && fluid_neighbor(point, cs) && match_level(point, cs);
}


/**
centroid_point returns the area of the interfrace fragment in a give cell. It
takes in the volume fraction field cs and fills mp with the interfacial 
centroid in the *global* coordinate system.

Note here n is the inward facing normal normalized so |n.x| + |n.y| + |n.z| = 1 */

double centroid_point (Point point, scalar s, face vector sf, coord * mp, coord * n, double * alpha)
{
    coord cellCenter = {x, y, z};
    *n = facet_normal (point, s, sf);
    *alpha = plane_alpha (s[], *n);
    double area = plane_area_center (*n, *alpha, mp);

    foreach_dimension()
        mp->x = cellCenter.x + mp->x*Delta;
    return area;
}


/**
closest_interface identifies the cell which contains the closest point to the 
interface from a ghost cell. It fills talpha, tfc, and bioff with the fragment's
intercept, cell-center coordinates, and offset indices, resp. */

coord closest_interface (Point point, vector mps, scalar cs, vector normals, scalar alphas,
                         double * talpha, coord * tfc, PointIBM * bioff, coord nu)
{
    coord temp_midPoint, temp_fluidCell = {0,0}, cc = {x, y, z};
    coord n;
    PointIBM ptemp = {0,0,0};

    bool found = false;
    double offmin = 1e6;

    foreach_neighbor() {
        if (cs[] > 0 && cs[] < 1) {

            /**
            The user can pass in a known nu. This is allowed when recursively
            searching for the closest point, though it is not used anymore. */

            n = !empty_coord(nu)? nu: (coord){normals.x[], normals.y[], normals.z[]};
            assert(!empty_coord(n));

            coord nplic = (coord){normals.x[], normals.y[], normals.z[]};
            double alpha = alphas[];

            /**
            Calculate the closest point to the current interface fragment. This
            is equivalent to finding the boundary intercept (BI) for a ghost cell. */

            coord poff = {(cc.x - x)/Delta, (cc.y - y)/Delta, (cc.z - z)/Delta}, bi = {0,0,0};
            double t = (alpha - dot_product(poff, nplic)) / (dot_product(n, nplic));
            foreach_dimension()
                bi.x = poff.x + t*n.x;

            coord big = {0, 0, 0}, fc = {x, y, z}; // fluid cell

            /**
            Put the boundary intercept in global coordinates. */
            foreach_dimension()
                big.x = fc.x + (bi.x*Delta);

            /**
            Make sure the boundary intercept lies within the cell. */
            if (fabs(bi.x) > 0.5 || fabs(bi.y) > 0.5 || fabs(bi.z) > 0.5)
                continue;

            /**
            We minimize the distance from the BI to the (ghost cell) cell center. */

            found = true;
            double offmin_temp = distance_coord(big, cc)/Delta;

            if (offmin_temp < offmin) {
                offmin = offmin_temp;
                temp_midPoint = (coord) {mps.x[], mps.y[], mps.z[]};
                temp_fluidCell = (coord){x, y, z};

                *talpha = alpha;
                ptemp = (PointIBM){(int)((x - cc.x)/Delta), (int)((y - cc.y)/Delta), (int)((z - cc.z)/Delta)};
            }
        }
    }

    /**
    If no boundary intercept is found, then the closest point to the interface is
    found using interfacial mid points.*/
    if (!found) {
      double min_distance = 1e6;

      foreach_neighbor() {
        double dx = mps.x[] - cc.x;
        double dy = mps.y[] - cc.y;
        double dz = mps.z[] - cc.z;

        n = !empty_coord(nu)? nu: (coord){normals.x[], normals.y[], normals.z[]};
       
        double alpha = alphas[];

        if ((mps.x[] || mps.y[] || mps.z[]) && 
             distance3D(dx, dy, dz) < min_distance) {

            temp_midPoint = (coord) {mps.x[], mps.y[], mps.z[]};

            temp_fluidCell = (coord){x, y, z};

            *talpha = alpha;

            min_distance = distance3D(dx, dy, dz);
            ptemp = (PointIBM){(int)((x - cc.x)/Delta), (int)((y - cc.y)/Delta), (int)((z - cc.z)/Delta)};
        }
      }
    }

    *tfc = temp_fluidCell;
    *bioff = ptemp;

    return temp_midPoint;
}


/**
# Normal interpolation

The normal vector at the boundary intercept mut be interpolated when it is not
readily available through analytical expressions. This function uses simple
linear interpolation to do so, though it is limited to 2D and is no longer used.
It is kept for documentation. */

coord interpolate_normal (Point point, coord bi, PointIBM bioff, vector mps, vector ns)
{
#if dimension == 2
    coord n1 = {ns.x[bioff.i, bioff.j, bioff.k], 
                ns.y[bioff.i, bioff.j, bioff.k],
                ns.z[bioff.i, bioff.j, bioff.k]};
    coord mp1 = {mps.x[bioff.i, bioff.j, bioff.k], 
                 mps.y[bioff.i, bioff.j, bioff.k],
                 mps.z[bioff.i, bioff.j, bioff.k]};
    
    PointIBM ioff = {sign(bi.x - mp1.x), sign(bi.y - mp1.y), sign(bi.z - mp1.z)};
    
    coord n2 = {ns.x[bioff.i + ioff.i, bioff.j, bioff.k], 
                ns.y[bioff.i + ioff.i, bioff.j, bioff.k],
                ns.z[bioff.i + ioff.i, bioff.j, bioff.k]};
    
    coord mp2 = {mps.x[bioff.i + ioff.i, bioff.j, bioff.k], 
                 mps.y[bioff.i + ioff.i, bioff.j, bioff.k],
                 mps.z[bioff.i + ioff.i, bioff.j, bioff.k]};

    coord n2t = {ns.x[bioff.i, bioff.j + ioff.j, bioff.k], 
                 ns.y[bioff.i, bioff.j + ioff.j, bioff.k],
                 ns.z[bioff.i, bioff.j + ioff.j, bioff.k]};
    
    coord mp2t ={mps.x[bioff.i, bioff.j + ioff.j, bioff.k], 
                 mps.y[bioff.i, bioff.j + ioff.j, bioff.k],
                 mps.z[bioff.i, bioff.j + ioff.j, bioff.k]};

    if (!empty_coord(n2) && !empty_coord(n2t)) {
        double d1 = distance_coord (mp2, bi);
        double d2 = distance_coord (mp2t, bi);
        if (d2 < d1) {
            mp2 = mp2t;
            n2 = n2t;
        }
    }
    else if (empty_coord(n2)) {
        mp2 = mp2t;
        n2 = n2t;
    }

    coord n = {0};
    if (empty_coord(n2) || equals_coord(n1, n2)) 
        n = n1;
     else {

        normalize(&n1);
        normalize(&n2);

        coord t1 = {-n1.y, n1.x};

        double s1 = dot_product(mp1, t1);
        double s2 = dot_product(mp2, t1);
        double sb = dot_product(bi, t1);

        bool check = (sb > s1 && sb < s2) || (sb > s2 && sb < s1) || (sb == s1) || (sb == s2);
        assert(check);

        double m = (sb - s1) / (s2 - s1 + SEPS);

        foreach_dimension()
            n.x = n1.x + m*(n2.x - n1.x);

        normalize_sum(&n); 
     }

     normalize_sum(&n);

     return n;

#else // dimension == 3

    coord n1 = {ns.x[bioff.i, bioff.j, bioff.k], 
                ns.y[bioff.i, bioff.j, bioff.k],
                ns.z[bioff.i, bioff.j, bioff.k]};

    return n1;

#endif // dimension == 3
}


/**
This function calculates the normal at a given pount (BI) using weighted
(or non-weighted) Least Squares (LSQ). 

Note, calculating s and r by doing $mp - bi$ makes s = r = 0 at the boundary 
intercept, which means we only need to solve for the a coefficient given the basis

n(s,r) = a + b*s + c*r

where s and r are the cartesian coordinates projected onto the tangent vectors 
t1 and t2, resp.

TODO: weight n by surface area of fragment. Larger surface area = less noise in n. */

coord interpolate_normal_lsq(Point point, coord bi, PointIBM bioff, vector mps,
                             vector ns)
{
    bool quad = false;
    bool weighted = true;

    coord n_anchor = {ns.x[bioff.i, bioff.j, bioff.k],
                      ns.y[bioff.i, bioff.j, bioff.k],
                      ns.z[bioff.i, bioff.j, bioff.k]};
    coord t1, t2;
    normal_and_tangents(&n_anchor, &t1, &t2);

    // Accumulators for LSQ fit n.x(s) = a_x + b_x*s and n.y(s) = a_y + b_y*s
    double Ss = 0, Sss = 0;
    double Snx = 0, Ssnx = 0, Sny = 0, Ssny = 0;
    double Kwi = 0;
    int K = 0;

#if dimension == 3
    double Sr = 0, Ssr = 0, Srr = 0;
    double Snz = 0, Ssnz = 0, Srnz = 0;
    double Srnx = 0, Srny = 0;
#endif

    coord gc = {x, y, z};

    double Ssss = 0, Ss4 = 0, Sssnx = 0, Sssny = 0;

    foreach_neighbor() {

        /**
        Exclude small interface fragments that make the interpolation noisy. */
        if (cs[] > 0 && cs[] < 1 - 1e-4) {

            coord mp_i = {mps.x[], mps.y[], mps.z[]};
            coord n_i  = {ns.x[],  ns.y[],  ns.z[]};
            normalize(&n_i);

            if (dot_product(n_i, n_anchor) < 0)
                foreach_dimension() n_i.x *= -1;

            double si = dot_product(subtract_coord(mp_i, bi), t1)/Delta;

            /**
            Weighted LSQ is used only if the function parameter weighted = true. 
            The weights are inversely proportional to the distance from the
            interfacial mid point and the ghost cell. */
            double wi = weighted? 1./(sq(distance_coord(mp_i, gc)/Delta) + sq(dimension)): 1;

            Ss  += si*wi;       Sss  += si*si*wi;
            Snx += n_i.x*wi;    Ssnx += si*n_i.x*wi;
            Sny += n_i.y*wi;    Ssny += si*n_i.y*wi;

        #if dimension == 3
            double ri = dot_product(subtract_coord(mp_i, bi), t2)/Delta;

            Sr   += ri*wi;        Ssr  += si*ri*wi;
            Srr  += ri*ri*wi;     Snz  += n_i.z*wi;
            Ssnz += si*n_i.z*wi;  Srnz += ri*n_i.z*wi;
            Srnx += ri*n_i.x*wi;  Srny += ri*n_i.y*wi;
        #endif

            if (quad) {
                Ssss  += si*si*si;
                Ss4   += si*si*si*si;
                Sssnx += si*si*n_i.x;
                Sssny += si*si*n_i.y;
            }

            Kwi += wi;
            K++;
        }
    }

    coord n;
    
    if (K < 2) {
        fprintf(stderr, "WARNING: not enough cells for LSQ normal interpolation, %d < 2\n", K);
        n = n_anchor;
    }
    else if (quad && K >= 3) {

        coord A1 = {Kwi, Ss, Sss};
        coord A2 = {Ss, Sss, Ssss};
        coord A3 = {Sss, Ssss, Ss4};
        double det = determinant_cols(A1, A2, A3);

        if (fabs(det) < 1e-10) {
            fprintf(stderr, "WARNING quadratic |det| = |%g| < %g\n", det, 1e-10);
            n = n_anchor;
        } 
        else {
            coord Bx = {Snx, Ssnx, Sssnx};
            double ax = determinant_cols(Bx, A2, A3) / det;

            coord By = {Sny, Ssny, Sssny};
            double ay = determinant_cols(By, A2, A3) / det;
            n = (coord){ax, ay, 0};
        }
    }
    else {
#if dimension == 2
        double det = Kwi*Sss - Ss*Ss;

        /**
        Tangent projection nearly collapses, fall back to anchor. */
        if (fabs(det) < 1e-10) {
            fprintf(stderr, "WARNING linear |det| = |%g| < 1e-10\n", det);
            n = n_anchor;
        } 
        else {
            double ax = (Sss*Snx - Ss*Ssnx) / det;
            double ay = (Sss*Sny - Ss*Ssny) / det;
            n = (coord) {ax, ay, 0};
       }
#else // dimension == 3
        coord ATA1 = {Kwi, Ss,  Sr};
        coord ATA2 = {Ss,  Sss, Ssr};
        coord ATA3 = {Sr,  Ssr, Srr};
        double det = determinant_cols(ATA1, ATA2, ATA3);

        if (fabs(det) < 1e-10) {
            fprintf(stderr, "WARNING quadratic |det| = |%g| < %g\n", det, 1e-10);
            n = n_anchor;
        } 
        else {
            coord Bx = {Snx, Ssnx, Srnx};
            double ax = determinant_cols(Bx, ATA2, ATA3) / det;

            coord By = {Sny, Ssny, Srny};
            double ay = determinant_cols(By, ATA2, ATA3) / det;

            coord Bz = {Snz, Ssnz, Srnz};
            double az = determinant_cols(Bz, ATA2, ATA3) / det;

            n = (coord){ax, ay, az};
        }
#endif // dimension == 3
    }

    normalize_sum(&n);
    return n;
}


/**
# Boundary intercept calculation

The function below returns the boundary intercept coordinate given a normal
vector, intercept, fluid cell center coordinates, and the volume fraction
field. 

Two normal vectors can be fined: ngc, the normal at the boundary intercept used
for ghost cells, and nplic, the normal of the PLIC reconstructed interface. One
can use this to find the boundary intercept on the PLIC interface using the 
analytical normal vector. */

coord boundary_int (Point point, coord ngc, coord nplic, double alpha, coord fc, scalar cs)
{
    coord poff = {(x - fc.x)/Delta, (y - fc.y)/Delta, (z - fc.z)/Delta}, bi = {0,0,0};
    double t = (alpha - dot_product(poff, nplic)) / (dot_product(ngc, nplic) + SEPS);
    foreach_dimension() {
        bi.x = poff.x + t*ngc.x;
        bi.x = fc.x + bi.x*Delta;
    }
    return bi;
}


coord direction_vector (coord p0, coord p1)
{
    return (coord){p1.x - p0.x, p1.y - p0.y, p1.z - p0.z};
}


/**
# Image point calculation

image_point takes in the coordinates of the boundary intercept and ghost cell and
returns the coordinates of the image point. */

coord image_point (coord bi, coord gc, coord n = {0})
{
     coord ip = {0};
     
     if (empty_coord(n)) {
        foreach_dimension()
            ip.x = gc.x + 2*(bi.x - gc.x);
     }
     else {
        double dist = distance_coord(bi, gc);
        normalize2(&n);

        foreach_dimension()
            ip.x = gc.x + 2*dist*n.x;
     }
     return ip;
}


/**
Similarly to image_point, fresh_image_point calculates the imagepoint for 
"fresh cells". This is only used for moving boundaries. */

coord fresh_image_point (coord boundaryInt, coord freshCell)
{
     double dx = freshCell.x - boundaryInt.x;
     double dy = freshCell.y - boundaryInt.y;
     double dz = freshCell.z - boundaryInt.z;

     coord imagePoint = {freshCell.x + dx, freshCell.y + dy, freshCell.z + dz};

     return imagePoint;
}


/**
# Interpolation at the image point

borders_boundary checks to see if a cell is adjacent to a domain boundary. 
Returns true if it has one or more neighbors that are inside the boundary.

Optionally, the user can provide two integer pointers (useri and userj) to be
filled with the indexs of the boundary neighbor. */

bool borders_boundary (Point point, int * useri = NULL, int * userj = NULL, int * userk = NULL)
{
#if TREE
    /**
    Look at directly adjacent neighbors (4 in 2D). */
    for (int d = 0; d < dimension; d++) {
	    for (int kk = -1; kk <= 1; kk += 2) {
            int _i = 0, _j = 0, _k = 0;
	        if (d == 0)
                _i = kk; 
            else if (d == 1)
                _j = kk; 
            else if (d == 2)
                _k = kk;
            /**
            check to see if neighboring cell is inside boundary. */
            if (neighbor(-_i,-_j,-_k).pid < 0) {

                if (useri) *useri = -_i;
                if (userj) *userj = -_j;
                if (userk) *userk = -_k;

                return true;
            }
            (void) _i; (void) _j; (void) _k; // to prevent unused variable warning
        }
    }
#endif
    return false;
}


#if _MPI
/**
borders_boundary checks to see if a cell is adjacent to a MPI boundary. 
Returns true if it has one or more neighbors that are inside the boundary. 

TODO: this can probably be simplified using foreach_direct_neighbor(). */

bool borders_mpi_boundary (Point point)
{
#if TREE 
    /**
    Look at directly adjacent neighbors (4 in 2D, 6 in 3D). */
    for (int d = 0; d < dimension; d++) {
	    for (int kk = -1; kk <= 1; kk += 2) {
            int i = 0, j = 0, k = 0;
	        if (d == 0)
                i = kk; 
            else if (d == 1)
                j = kk; 
            else if (d == 2)
                k = kk;
            /**
            check to see if neighboring cell is outside the current MPI 
            rank/processor. */
            if (neighbor(-i,-j,-k).pid != pid())
                return true;
            (void) k;
        }
    }
#endif
    return false;
}
#endif


/**
image_offsets fills integers xOffset and yOffset with the index of the cell containing
the given image point w.r.t the ghost cell's stencil.

e.g., if an image point is 2 cells to the left of it's respective ghost cell, then
xOffset = -2 and yOffset = 0. */

int image_offsets (Point point, coord imagePoint, int *xOffset, int *yOffset, 
                   int *zOffset = NULL)
{
    coord ghostCell = {x,y,z};
    int offset_x = 0, offset_y = 0, offset_z = 0;
    foreach_dimension() {
        double d = fabs(imagePoint.x - ghostCell.x) / Delta;
        int dsign = sign (imagePoint.x - ghostCell.x);
        if (d >= 1.5) {
            if (d >= 2) 
              offset_x = 3 * dsign;
            else
              offset_x = 2 * dsign;
        }
        else if (d >= 0.5) {
            offset_x = 1 * dsign;
        }
    }
    if (xOffset != NULL && yOffset != NULL) {
        *xOffset = offset_x;
        *yOffset = offset_y;
    }
    if (zOffset != NULL) {
        *zOffset = offset_z;
    }

    return 1;
}


/**
This function fills pints with the coordinates of the bounding nodes used for
interpolation and pnodes with their respective indices within the stencil. */

void get_interpolation_points (Point point, const int m, coord pints[m], 
                               PointIBM pnodes[m], PointIBM poff, PointIBM pnode)
{
    coord icell = {x + Delta*poff.i, y + Delta*poff.j, z + Delta*poff.k};

    pints[0]  = (coord){icell.x,                 icell.y,                 icell.z};
    pints[1]  = (coord){icell.x + pnode.i*Delta, icell.y,                 icell.z};
    pints[2]  = (coord){icell.x + pnode.i*Delta, icell.y + pnode.j*Delta, icell.z};
    pints[3]  = (coord){icell.x,                 icell.y + pnode.j*Delta, icell.z};

    pnodes[0] = (PointIBM){0,       0,       0};
    pnodes[1] = (PointIBM){pnode.i, 0,       0};
    pnodes[2] = (PointIBM){pnode.i, pnode.j, 0};
    pnodes[3] = (PointIBM){0,       pnode.j, 0};

#if dimension == 3
    pints[4]  = (coord){icell.x,                 icell.y,                 icell.z + pnode.k*Delta};
    pints[5]  = (coord){icell.x + pnode.i*Delta, icell.y,                 icell.z + pnode.k*Delta};
    pints[6]  = (coord){icell.x + pnode.i*Delta, icell.y + pnode.j*Delta, icell.z + pnode.k*Delta};
    pints[7]  = (coord){icell.x,                 icell.y + pnode.j*Delta, icell.z + pnode.k*Delta};

    pnodes[4] = (PointIBM){0,       0,       pnode.k};
    pnodes[5] = (PointIBM){pnode.i, 0,       pnode.k};
    pnodes[6] = (PointIBM){pnode.i, pnode.j, pnode.k};
    pnodes[7] = (PointIBM){0,       pnode.j, pnode.k};
#endif

}

#define rows (1 << dimension)
#define cols (rows + 1)


/**
This function fills rmatrix with the row used in the interpolation system
corrseponding to a entry in the augmented Vandermonde matrix.

It considers boundary conditions when a bounding node is a ghost cell to close
the system. */

void fluid_only (Point point, scalar s, const int n, double rmatrix[n],
                  PointIBM poff, PointIBM pnode, PointIBM pbound, 
                  coord * pcell, double sval, coord ipoint,
                  vector mps, vector ns, scalar alphas)
{
    double val = sval;

    int xx = poff.i + pnode.i, yy = poff.j + pnode.j, zz = poff.k + pnode.k;
    int type = 0;

    /**
    If a node is a ghost cell, then the cell's entry in the linear system is
    changed based on the wall's boundary condition. */

    if (cs[xx,yy,zz] <= GCV && cs[xx,yy,zz] > 0.) {
        *pcell = (coord){mps.x[xx,yy,zz], mps.y[xx,yy,zz], mps.z[xx,yy,zz]};

        /**
        Move point more if cell is inside domain boundary. */

        if (xx == pbound.i) pcell->x += pbound.i*Delta;
        if (yy == pbound.j) pcell->y += pbound.j*Delta;
        if (zz == pbound.k) pcell->z += pbound.k*Delta;

        /**
        Check the boundary condition. */
        bool bctype[2] = {false, false};
        double bc = s.boundary[immersed] (point, point, s, bctype);

        bool dirichlet = bctype[0], navierslip = bctype[1];

        if (dirichlet && !navierslip) {
            val = bc;
            type = 0;
        }
        else if (dirichlet && navierslip) {
            val = bc; // assumes stationary wall
            type = 2;
        }
        else { // neumann
            val = bc;
            type = 1;
        }
    }

    /**
    The row in the linear system is now constructed based on the coordinates and
    boundary condition (if the node is a ghost cell). */

    if (type == 0) { // normal or dirichlet
#if dimension == 2
        
        memcpy(rmatrix, (double[]){pcell->x*pcell->y, pcell->x, pcell->y, 1, val}, cols*sizeof(double));
#else // dimension == 3
        memcpy(rmatrix, (double[]){pcell->x*pcell->y*pcell->z,
                                   pcell->x*pcell->y, pcell->x*pcell->z, pcell->y*pcell->z,
                                   pcell->x, pcell->y, pcell->z, 1, val}, cols*sizeof(double));
#endif // dimension == 2
    }
    else if (type == 1) { // neumann
        coord n = {ns.x[xx,yy,zz], ns.y[xx,yy,zz], ns.z[xx,yy,zz]};
        normalize(&n);
#if dimension == 2
        memcpy(rmatrix, (double[]){n.x*pcell->y + n.y*pcell->x, n.x, n.y, 0, val}, cols*sizeof(double));
#else // dimension == 3
        memcpy(rmatrix, (double[]){n.x*pcell->y*pcell->z + n.y*pcell->x*pcell->z + n.z*pcell->x*pcell->z,
                                   n.x*pcell->y + n.y*pcell->x,
                                   n.x*pcell->z + n.z*pcell->x,
                                   n.y*pcell->z + n.z*pcell->y,
                                   n.x, n.y, n.z, 0, val}, cols*sizeof(double));
#endif // dimension == 2
    }
    else if (type == 2) { // navier-slip
        coord n = {ns.x[xx,yy,zz], ns.y[xx,yy,zz], ns.z[xx,yy,zz]};
        double mag = magnitude_coord(n) + SEPS;
        normalize(&n);
        double alpha = alphas[xx,yy,zz];

        coord gc = {x + xx*Delta, y + yy*Delta, z + zz*Delta}, bi; // ghost cell, boundary intercept
        foreach_dimension()
            bi.x = alpha/mag * n.x;

        foreach_dimension()
            pcell->x = gc.x + bi.x*Delta;

        coord d = direction_vector(gc, *pcell);
        double dmag = magnitude_coord(d);        // distance from ghost cell to boundary intercept
        double term = (1 - dmag/(val + SEPS));
        double usolid = 0.0*dmag/(val + SEPS);  // navier b.c only works for stationary solids right now

#if dimension == 2
        memcpy(rmatrix, (double[]){gc.x*gc.y - (pcell->x*pcell->y)*term,
                                   gc.x - pcell->x*term,
                                   gc.y - pcell->y*term,
                                   1 - term, usolid}, cols*sizeof(double));
#else // dimension == 3
        memcpy(rmatrix, (double[]){gc.x*gc.y*gc.z - pcell->x*pcell->y*pcell->z*term,
                                   gc.x*gc.y - pcell->x*pcell->y*term,
                                   gc.x*gc.z - pcell->x*pcell->z*term,
                                   gc.y*gc.z - pcell->y*pcell->z*term,
                                   gc.x - pcell->x*term,
                                   gc.y - pcell->y*term,
                                   gc.z - pcell->z*term,
                                   1 - term, usolid}, cols*sizeof(double));
#endif // dimension == 2
    }

    NOT_UNUSED(zz); // to prevent unused variable warning
}


/**
This function constructs the augmented Vandermonde matrix, using the fluid_only
function to construct each row. */
void get_interpolation_matrix (Point point, int m, int n, double matrix[m][n],
                               double snodes[m], scalar s, PointIBM poff, 
                               PointIBM pnode, PointIBM pbound, coord ip, 
                               vector mps, vector ns, scalar alphas)
{
    /**
    Calculate (global) coordinates and relative indices of interpolation cells. */
    coord pints[rows];
    PointIBM pnodes[rows];
    get_interpolation_points(point, m, pints, pnodes, poff, pnode);

    /**
    Assemble matrix row by row. */
    for (int row = 0; row < m; ++row) {

        /**
        If a cell for interpolating is a ghost cell, move the point to the
        interface and change the velocity to the correct boundary condition. */
        fluid_only(point, s, n, matrix[row], poff, pnodes[row], pbound, 
                   &pints[row], snodes[row], ip, mps, ns, alphas);
    }
}

/**
Fill nodes with the bounding values of s to be used in the interpolation system.*/
void fill_nodal_array (Point point, scalar s, double nodes[rows], 
                       PointIBM off, PointIBM d)
{
    nodes[0] = s[off.i,       off.j,       off.k];
    nodes[1] = s[off.i + d.i, off.j,       off.k];
    nodes[2] = s[off.i + d.i, off.j + d.j, off.k];
    nodes[3] = s[off.i,       off.j + d.j, off.k];
#if dimension == 3
    nodes[4] = s[off.i,       off.j,       off.k + d.k];
    nodes[5] = s[off.i + d.i, off.j,       off.k + d.k];
    nodes[6] = s[off.i + d.i, off.j + d.j, off.k + d.k];
    nodes[7] = s[off.i,       off.j + d.j, off.k + d.k];
#endif
}

/**
This function constructs and solves the interpolation matrix and returns the
value of s at the image point, ip.*/
double solve_interpolation_ibm (Point point, scalar s, double nodes[rows],
                                PointIBM off, PointIBM d, PointIBM boff,
                                vector mps, vector ns, scalar alphas,
                                coord ip)
{
    /**
    Fill the interpolation matrix to be inverted. Here, the matrix is an
    agumented matrix of the Vandermonde matrix and the vector containg the values
    of s at the bounding nodes, which is snodes in this case. */

    double matrix[rows][cols];
    get_interpolation_matrix(point, rows, cols, matrix, nodes, s, 
                             off, d, boff, ip, mps, ns, alphas);

    /**
    Solve the linear system to get the coefficients */
    double coeff[rows];
    gauss_elim (rows, cols, matrix, coeff);

    /**
    Use the coefficients and the coordinate of the image point to find the
    interpolated value of s at it. */

#if dimension == 2
    double sip = coeff[0] * ip.x * ip.y +
                 coeff[1] * ip.x +
                 coeff[2] * ip.y +
                 coeff[3];
#else
    double sip = coeff[0] * ip.x * ip.y * ip.z +
                 coeff[1] * ip.x * ip.y +
                 coeff[2] * ip.x * ip.z +
                 coeff[3] * ip.y * ip.z +
                 coeff[4] * ip.x +
                 coeff[5] * ip.y +
                 coeff[6] * ip.z +
                 coeff[7];
#endif
    return sip;
}


/**
## Interpolation of a scalar field

This function handles the interpolation of a scalar field s and a given point
ip using bilinear (2D) or trilinear (3D) interpolation. */

double interpolate_scalar_ibm (Point point, scalar s, coord ip,
                              vector mps, vector ns, scalar alphas)
{
    /**
    If the ghost cell is near a boundary, we slightly adjust the interpolation
    scheme. */

    int boffx = 0, boffy = 0, boffz = 0; // boundary offsets
    borders_boundary (point, &boffx, &boffy, &boffz);
    PointIBM boff = {boffx, boffy, boffz};

    /**
    We calculate the coordinates of the cell of which image point lies in and 
    the location of it within the stencil of the current point, the later
    stored in the three offset variables below. */

    int xx = 0, yy = 0, zz = 0; 
    image_offsets (point, ip, &xx, &yy, &zz);
     PointIBM off = {xx, yy, zz};   

    assert (abs(xx) <= 2 && abs(zz) <= 2 && abs(zz) <= 2);

    coord ic = {x + Delta*xx, y + Delta*yy, z + Delta*zz};
    
    /**
    Calculate indices for interpolating stencil. */
    PointIBM d = {sign(ip.x - ic.x), sign(ip.y - ic.y), sign(ip.z - ic.z)};

    /**
    Fill the array storying the bounding points. */
    double nodes[rows];
    fill_nodal_array(point, s, nodes, off, d);

    /**
    Construct and solve the linear system to interpolate s at the image point. */
    return solve_interpolation_ibm (point, s, nodes, off, d, boff, mps, ns, alphas, ip);
}

/**
## Interpolation of a vector field

This function handles the interpolation of a vector field u and a given point
ip using bilinear (2D) or trilinear (3D) interpolation. */

coord interpolate_vector_ibm (Point point, vector u, coord ip, coord n, 
                              vector mps, vector ns, scalar alphas)
{
    /**
    If the ghost cell is near a boundary, we slightly adjust the interpolation
    scheme. */

    int boffx = 0, boffy = 0, boffz = 0; // boundary offsets
    borders_boundary (point, &boffx, &boffy, &boffz);
    PointIBM boff = {boffx, boffy, boffz};

    /**
    We calculate the coordinates of the cell of which image point lies in and 
    the location of it within the stencil of the current point, the later
    stored in the three offset variables below. */

    int xx = 0, yy = 0, zz = 0; 
    image_offsets (point, ip, &xx, &yy, &zz);
     PointIBM off = {xx, yy, zz};   

    // assert (abs(xx) <= 2 && abs(zz) <= 2 && abs(zz) <= 2);
    xx = clamp(xx, -2, 2);
    yy = clamp(yy, -2, 2);
    zz = clamp(zz, -2, 2);

    coord ic = {x + Delta*xx, y + Delta*yy, z + Delta*zz};
    
    /**
    Calculate indices for interpolating stencil. */
    PointIBM d = {sign(ip.x - ic.x), sign(ip.y - ic.y), sign(ip.z - ic.z)};

    /**
    Fill the arrays storying the bounding points. */
    double nodes_x[rows];
    fill_nodal_array(point, u.x, nodes_x, off, d);

    double nodes_y[rows];
    fill_nodal_array(point, u.y, nodes_y, off, d);

#if dimension == 3
    double nodes_z[rows];
    fill_nodal_array(point, u.z, nodes_z, off, d);
#endif

    /**
    Project the vectors to the local wall coordinate system: x->n, y->t, z->r. */
    coord t1, t2;
    normal_and_tangents (&n, &t1, &t2);

    for (int i = 0; i < rows; ++i) {
        coord ui = {0};
        foreach_dimension()
            ui.x = nodes_x[i];
        coord uiproj = {dot_product(ui, n),
                        dot_product(ui, t1),
                        dot_product(ui, t2)};
        foreach_dimension()
            nodes_x[i] = uiproj.x;
    }

    /**
    Construct and solve the linear system to interpolate s at the image point. */   
    coord uip = {0}; 
    uip.x = solve_interpolation_ibm(point, u.x, nodes_x, off, d, boff, mps, ns, alphas, ip);
    uip.y = solve_interpolation_ibm(point, u.y, nodes_y, off, d, boff, mps, ns, alphas, ip);   
#if dimension == 3 
    uip.z = solve_interpolation_ibm(point, u.z, nodes_z, off, d, boff, mps, ns, alphas, ip);   
#endif

    /**
    Project the interpolated vector value back to cartesian coordinates. */
    double gcn = uip.x, gct1 = uip.y, gct2 = uip.z;
    foreach_dimension()
        uip.x = gcn*n.x + gct1*t1.x + gct2*t2.x;

    return uip;
}

#undef rows
#undef cols


/**
# Support functions

This function is adapted from embed.h. It removes any cell with inconsistent 
volume/surface fractions. Additionally, it removes an cell with a fluid volume
fractions less than smin. */

struct Cleanup {
  scalar c;
  face vector s;
  double smin;   // minimum surface fraction (optional)
  bool opposite; // whether to eliminate 'thin tubes' (optional)
};

trace
int fractions_cleanup (scalar c, face vector s,
		       double smin = 0., bool opposite = false)
{
  
  /*
  Since both surface and volume fractions are altered, iterations are
  needed. This reflects the fact that changes are coupled spatially
  through the topology of the domain: for examples, long, unresolved
  "filaments" may need many iterations to be fully removed. */
  
  int changed = 1, schanged = 0, i;
  for (i = 0; i < 100 && changed; i++) {

    /**
    Face fractions of empty cells must be zero. */
   
    foreach_face() {
      if (s.x[] && ((!c[] || !c[-1]) || s.x[] < smin || 
                   ((c[] > 0 && c[] < smin && level == depth()) || 
                   (c[-1] > 0 && c[-1] < smin && level == depth())))) {
        s.x[] = 0.;
      }
    }

    changed = 0;
    foreach(reduction(+:changed))
      if (c[] > 0. && c[] < 1.) {
        int n = 0;
        foreach_dimension() {
          for (int i = 0; i <= 1; i++)
	    if (s.x[i] > 0.)
	      n++;

	  /**
	  If opposite surface fractions are zero (and volume fraction
	  is non-zero), then we are dealing with a thin "tube", which
	  we just remove because it can sometimes lead to
	  non-convergence when
	  [projecting](navier-stokes/centered.h#approximate-projection)
	  the velocity field. */

	  if (opposite && s.x[] == 0. && s.x[1] == 0.)
	    c[] = 0., changed++;
	}

	/**
	The number of "non-empty" faces (i.e. faces which have a
	surface fraction larger than epsilon) cannot be smaller than
	the dimension (the limiting cases correspond to a triangle in
	2D and a tetrahedron in 3D). */
	
	if (n < dimension)
	  c[] = 0., changed++;
      }

    schanged += changed;
  }
  if (changed)
    fprintf (stderr, "WARNING: fractions_cleanup() did not converge after "
	     "%d iterations\n", i);
  return schanged;
}


/**
ibm_geometry is taken from embed's embed_geometry. Given a point containing a
interface fragment, the function returns the area of the fragment and fills p and
n with the coordinates of the interfacial mid point and normal vector, respectively.

Note, n is normalized differently from the MYC's way. Using the interface_normal
function from fractions.h, |n.x| + |n.y| = 1. Here, n is normalized using the vectors
magnitude, i.e. sqrt( sq(n.x) + sq(n.y) ) = 1.

Note area must be multiplied by the cell length (or area in 3D) to get in terms of
phyiscal "units". Also, the midpoint, p, is in the cell's local coordinate system. */

static inline
double ibm_geometry (Point point, coord * p, coord * n)
{
    *n = facet_normal (point, cs, fs);
    double alpha = plane_alpha (cs[], *n);

    double area = plane_area_center (*n, alpha, p);
    foreach_dimension()
        n->x *= -1;
    normalize (n);

    return area;
}


static inline
double ibm0_geometry (Point point, coord * p, coord * n, scalar cs1, face vector fs1)
{
    *n = facet_normal (point, cs1, fs1);
    double alpha = plane_alpha (cs1[], *n);
    double area = plane_area_center (*n, alpha, p);
    foreach_dimension()
        n->x *= -1;
    normalize (n);

    return area;
}


/**
The next few functions are taken from embed to calculate interfacial force.
*/

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
					   coord n, coord p, double bc,
					   double * coef)
{
  //foreach_dimension()
  //  n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }

  //fprintf(stderr, "(%g,%g) d0=%g d1=%g bc=%g v0=%g v1=%g cs=%g\n", x, y, d[0], d[1], bc, v[0], v[1], cs[]);

  if (v[0] == nodata) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient (Point point, scalar s, scalar cs,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, cs, n, p, bc, coef);
#endif // dimension == 3
  return nodata;
}

static inline
coord ibm_gradient (Point point, vector u, coord p, coord n)
{
    coord dudn;
    foreach_dimension() {
        bool dirichlet = false;
        double vb = u.x.boundary[immersed] (point, point, u.x, &dirichlet);
        if (dirichlet) {
            double val;
            dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, &val);
        }
        else
            dudn.x = vb;
        if (dudn.x == nodata)
          dudn.x = 0.;
    }
    return dudn;
}

double ibm_vorticity (Point point, vector u, coord p, coord n)
{
    coord dudn = ibm_gradient (point, u, p, n);

    return -(dudn.y*n.x - dudn.x*n.y);
}


/**
ibm_force calculates the pressure force, Fp, and viscous force, Fmu, acting on the
immersed boundary. */

void ibm_force (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};
    foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        // if cell contains boundary intercept
        if (cs[] > 0. && cs[] < 1.) {
            coord midPoint, n, b;
            double area = ibm_geometry (point, &b, &n);
#if AXI
            double val = cm[];
#else
            double val = 1;
#endif
            area *= val*pow (Delta, dimension - 1); // is cm[]* right.. for axi?

            
            coord cellCenter = {x,y,z};
            foreach_dimension() {
                midPoint.x = cellCenter.x + b.x*Delta;
            }
            // calculate pressure force
            double boundaryPressure = extrapolate_scalar (point, cs, midPoint, n, p);
            double Fn = area * boundaryPressure;

            foreach_dimension()
                Fps.x -= Fn * n.x;

            // calculate shear force
            if (constant(mu.x) != 0.) {
            	double mua = 0., fa = 0.;

            	foreach_dimension() {
                    mua += mu.x[] + mu.x[1];
                    fa  += fm.x[] + fm.x[1];
	            }

                mua /= (fa + SEPS);
                coord velocityGrad = ibm_gradient (point, u, b, n);

#if dimension == 2
                foreach_dimension()
                    Fmus.x -= area * mua* (velocityGrad.x * (sq(n.x) + 1.) + 
                                           velocityGrad.y * -n.x * -n.y);
#else
                foreach_dimension()
                    Fmus.x -= area * mua * (velocityGrad.x * (sq(n.x) + 1.) + 
                                            velocityGrad.y * -n.x * -n.y +
                                            velocityGrad.z * -n.x * -n.z);
#endif
            }
        }
    }
    *Fp = Fps;
    *Fmu = Fmus;
}


/**
This function fills cf with the skin friction coefficient for each interfacial cell.

TODO: this assume no-slip condition; make it work for all B.C! */

double skin_friction (vector u, face vector mu, scalar cf)
{
    double cftotal = 0;
    foreach (reduction(+:cftotal)) {
        if (cs[] > 0 && cs[] < 1) {
            coord n, b;
            ibm_geometry (point, &b, &n);
            
            // calculate shear force
            double mua = 0., fa = 0.;

            foreach_dimension() {
                mua += mu.x[] + mu.x[1];
                fa  += fm.x[] + fm.x[1];
	        }

            mua /= (fa + SEPS);
            coord dudn = ibm_gradient (point, u, b, n);
            coord tau = {0,0,0};
           
#if dimension == 2
            foreach_dimension()
                tau.x -= mua* (dudn.x * (sq(n.x) + 1.) + 
                               dudn.y * -n.x * -n.y);
            cf[] = distance(tau.x, tau.y);
#else
            foreach_dimension()
                tau.x -= mua * (dudn.x * (sq(n.x) + 1.) + 
                                dudn.y * -n.x * -n.y +
                                dudn.z * -n.x * -n.z);
            cf[] = magnitude_coord(tau);
#endif
            cftotal += cf[];
        }
        else
            cf[] = 0;
    }
    return cftotal;
}


/**
This function is used in the BCG advection scheme in my-bcg.h. */
static void gradients_ibm (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
      if (s.gradient)
	foreach_dimension() {
      if (!gc[])
        v.x[] = 0.;
      else
	    v.x[] = s.gradient (s[-1], s[], s[1])/Delta;
	}
      else // centered
	    foreach_dimension() {
          if (!gc[])
            v.x[] = 0.;
          else
	        v.x[] = (s[1] - s[-1])/(2.*Delta);
	    }
    }
  }
}


/**
# Face centered interpolation 

The face_value macro is overridden to account for the immersed boundary. It is
only used when calculating the face velocity field uf. */

#define ibm_avg(a,i,j,k)							\
  ((a[i,j,k]*(0.1 + cs[i,j,k]) + a[i-1,j,k]*(0.1 + cs[i-1,j,k]))/	\
   (cs[i,j,k] + cs[i-1,j,k] + 0.2))

#undef face_value
#define face_value(a,i)							\
  (fs.x[i] < 1. && fs.x[i] > 0. ?				\
   0.5*(a[i] + a[i-1]):                         \
   ibm_avg(a,i,0,0))

/**
# Metric event

The metric event is used to set the metric fields, fm and cm, to the gcf and
gc fields, respectively. It is also used to specifiy the prolongation and
refinement operations for each field, which are defined in ibm-tree.h. */

#if TREE
#include "ibm-tree.h"
#endif

event metric (i = 0)
{
    if (is_constant (fm.x)) {
        foreach_dimension()
            assert (constant (fm.x) == 1.);
        fm = gcf;
    }
    foreach_face() {
        gcf.x[] = 1.;
        fs.x[] = 1;
    }
    if (is_constant (cm)) {
        assert (constant (cm) == 1.);
        cm = gc;
    }
    foreach() 
        gc[] = cs[] = csp[] = 1;

#if TREE
    /**
    Set prolongation and refining functions. */

    gc.refine = fraction_refine_metric;
    gc.prolongation = fraction_refine_metric;

    cs0.refine = cs.refine = ibm_fraction_refine;
    cs0.prolongation = cs.prolongation = fraction_refine;

    csp.refine = csp.prolongation = fraction_refine;

    foreach_dimension() {
        gcf.x.prolongation = refine_metric_injection_x;
        fs.x.prolongation = ibm_face_fraction_refine_x;
        gcf.x.restriction = restriction_face_metric; // is this really necessary?
    }
#endif
    restriction ({cs, csp, fs, gcf, gc});

    boundary(all);
}

