#define dv2() (pow(Delta,dimension)*ibm[]*cm[])

#include "fractions.h"

#define BGHOSTS 2
#define IBM 1
#define LIMIT 1e100
#define INT_TOL 1e-7    // tolerance used for volume fraction fields (interface tolerance)
#define VTOL 1e-11

#undef SEPS
#define SEPS 1e-30

scalar ibm[];
scalar ibm0[];          // solid volume fraction field of previous timestep
face vector ibmf[];
face vector ibmf0[];

// metric fields
scalar ibmCells[];
face vector ibmFaces[];

double (* metric_ibm_factor) (Point, coord) = NULL; // for axi

/**
when true, immersed BCs will respect the n and t (and r) notation, otherwise,
they will be decayed to n ≡ x, t ≡ y, and r ≡ z 
*/
bool local_bc_coordinates = true; 

#if 0
typedef struct solidVelo {
    void (* move_solid_x) (scalar * ibm, face vector * ibmf);
    void (* move_solid_y) (scalar * ibm, face vector * ibmf);
} solidVelo;

solidVelo usolid;
#endif

typedef struct fragment {
    coord n;
    double alpha;
    double c;  // solid volume fraction field (ibm)
} fragment;


void fill_fragment (double c, coord n, fragment * frag)
{
    frag->c = c;
    frag->n = n;
    frag->alpha = plane_alpha (c, n);
}

typedef struct PointIBM {
    int i, j, k;
} PointIBM;


#define distance(a,b) sqrt(sq(a) + sq(b))
#define distance3D(a,b,c) sqrt(sq(a) + sq(b) + sq(c))

#define on_interface(a) (a[] > 0+INT_TOL && a[] < 1-INT_TOL)
//#define on_interface(a,_TOL) (a[] > 0+_TOL && a[] < 1.-_TOL)

#define is_mostly_solid(a, i) (a[i] > 0+INT_TOL && a[i] <= 0.5)
#define is_fresh_cell(a0, a) (a0[] <= 0.5 && a[] > 0.5)

bid immersed;

attribute {
    vector mp; // boundary intercepts (TODO: should name this better)
}

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
  return data ? ((bool *)data)[0] = true, 0 : - s[];
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
  return data ? ((bool *)data)[0] = false, 0 : s[];
}

macro2
double navier_slip (double expr, Point point = point,
		  scalar s = _s, bool * data = data)
{
  return data ? ibm_area_center (point, s, &x, &y, &z),
    ((bool *)data)[0] = true, ((bool *)data)[1] = true, expr : 2.*expr - s[];
}

coord cross_product (coord a, coord b)
{
    coord c = {
        (a.y*b.z) - (a.z*b.y) ,
      -((a.x*b.z) - (b.x*a.z)),
        (a.x*b.y) - (a.y*b.x)
    };
    return c;
}

double determinant (coord a, coord b)
{
    coord c = cross_product(a,b);
    return distance3D(c.x,c.y,c.z); // area
}


int approx_equal (coord p1, coord p2, double TOL = VTOL)
{
    return fabs(p1.x - p2.x) <= TOL && fabs(p1.y - p2.y) <= TOL && fabs(p1.z - p2.z) <= TOL;
}

int approx_equal_double (double a, double b, double TOL = VTOL)
{
    return fabs(a - b) <= TOL;
}


double dot_product_norm (coord a, coord b)
{
    normalize(&a); normalize(&b);

    double product = 0;
    foreach_dimension()
        product += a.x*b.x;

    return product;
}

double dot_product (coord a, coord b)
{
    double product = 0;
    foreach_dimension()
        product += a.x*b.x;

    return product;
}

double dot_product_angle (coord a, coord b)
{
    double product = clamp(dot_product(a, b), -1, 1);
    assert(fabs(product) <= 1);
    return acos(dot_product(a, b));
}

// normalize() but with SEPS in the denominator
void normalize2 (coord * n)
{
    double norm = 0;
    foreach_dimension()
        norm += sq(n->x);
    norm = sqrt(norm);
    foreach_dimension()
        n->x /= norm + SEPS;
}

void normalize_sum (coord * n)
{
    double norm = 0;
    foreach_dimension()
        norm += fabs(n->x);
    foreach_dimension()
        n->x /= norm + SEPS;
}

/*
normal_and_tangents accepts a normal vector and fill t1 and t2 with the 
corresponding tangent vector(s) (one in 2D and two in 3D), all normalized.
*/
void normal_and_tangents (coord * n, coord * t1, coord * t2)
{
    normalize2(n);
#if dimension == 2

    coord t1_tmp = {-n->y, n->x};
    *t1 = t1_tmp;
    *t2 = (coord){0,0,0};

#else // dimension == 3

    coord a = {0,0,1};
    if ((!fabs(dot_product(*n, a))) < 0.9) {
        a = (coord){1,0,0};
    }
    coord t1_tmp = cross_product(*n, a);
    double det = distance3D(t1_tmp.x, t1_tmp.y, t1_tmp.z); // determinant

    assert(fabs(det) > 1e-15);

    foreach_dimension()
        t1->x = t1_tmp.x/det;
    
    *t2 = cross_product(*n, *t1);

#endif // dimension == 3
}

/*
This function takes returns true if the given point has a direct neighbor that
has no liquid volume fraction, i.e. ibm == 0, and fills pc and n with the midpoint
and corresponding normal, respectively.

TODO: should only check neighbors sharing a face, N, S, E, or W.
*/

bool empty_neighbor (Point point, coord * pc, coord * n, scalar ibm)
{
    coord pc_temp, cellCenter = {x, y, z};
    double ibm_temp = ibm[];
    double max_d = 1e6;
    int neighbor = 0;

    foreach_neighbor(1) {
        double distance2Cell = distance3D(x - cellCenter.x, y - cellCenter.y, z - cellCenter.z);
        if (ibm[] == 0 && ibm_temp == 1 && distance2Cell < max_d) {
            pc_temp.x = (cellCenter.x + x) / 2.;
            pc_temp.y = (cellCenter.y + y) / 2.;
            pc_temp.z = (cellCenter.z + z) / 2.;

            max_d = distance3D(x - cellCenter.x, y - cellCenter.y, z - cellCenter.z);

            neighbor = 1;
            *pc = pc_temp;
        }
    }
   
    // Calculate the normal. n can only be {1,0}, {-1,0}, {0,1}, or {0,-1}.
    // should probably not compare floating point values like this.
    coord n_temp = {pc_temp.x != cellCenter.x, 
                    pc_temp.y != cellCenter.y,
                    pc_temp.z != cellCenter.z};
    foreach_dimension() {
        n_temp.x *= sign2(pc_temp.x - cellCenter.x);
    }

    *n = n_temp;

    return neighbor;
}


/*
Checks to see if the given point has at least 1 fluid neigbor touching one of
it's faces. If so, returns true.

TODO: change algorithm to only check neighbors, not the cell itself (ibm[0,0])
      - should only require using i += 2 instead of i++
*/

bool fluid_neighbor (Point point, scalar ibm)
{
    // check left and right neighbors
    for(int i = -1; i <= 1; i++)
        if (ibm[i] > 0.5)
            return true;

    // check top and bottom neighbors
    for(int j = -1; j <= 1; j++)
        if (ibm[0, j] > 0.5)
            return true;

#if dimension == 3
    // check front and back neighbors
    for(int k = -1; k <= 1; k++)
        if (ibm[0, 0, k] > 0.5)
            return true; 
#endif

    return false;
}


/*
match_level is used to make sure that the neighboring fluid cell (assuming there 
is one) does not contain children that are ghost cells which can undesireably lead 
to two layers of ghost cells and constant refining/coarsening.

TODO: only check N, S, E, and W neighbors, not entire 3x3 stencil
*/

bool match_level (Point point, scalar ibm)
{
    foreach_neighbor(1) {
        if (ibm[] > 0.5 && is_leaf(cell) && is_active(cell))
            return true;
    }
    return false;
}


/*
is_ghost_cell returns true if the given cell shares a face with a fluid cell,
ibm > 0.5, and the volume fraction is less than or equal to 0.5.
*/

bool is_ghost_cell (Point point, scalar ibm)
{
   return ibm[] <= 0.5 && fluid_neighbor(point, ibm) && match_level(point, ibm);
}


/*
centroid_point returns the area of the interfrace fragment in a give cell. It
takes in the volume fraction field ibm and fills midPoint with the interfacial 
centroid in the GLOBAL coordinate system.

Note here n is the inward facing normal normalized so |n.x| + |n.y| + |n.z| = 1
*/

double centroid_point (Point point, scalar ibm, coord * midPoint, coord * n, double * alpha)
{
    coord cellCenter = {x, y, z};
    *n = facet_normal (point, ibm, ibmf);
    *alpha = plane_alpha (ibm[], *n);
    double area = plane_area_center (*n, *alpha, midPoint);

    foreach_dimension()
        midPoint->x = cellCenter.x + midPoint->x*Delta;
    return area;
}


/**
reconstruction_ibm accepts an additional face vector field that represents the
face solid volume fraction (ibmf) to be used when calculating the normal, n.

TODO: what if interface perfectly cuts cell face? use interfacial() instead? must
      be as cheap as possible.
*/

void reconstruction_ibm (const scalar c, const face vector cf, vector n, scalar alpha)
{
    foreach() {
        if (c[] <= 0. || c[] >= 1.) {
            alpha[] = 0.;
            foreach_dimension()
                n.x[] = 0.;
        }
        else {
            coord m = facet_normal (point, c, cf);
            foreach_dimension()
                n.x[] = m.x;
            alpha[] = plane_alpha(c[], m);
        }
    }

#if TREE
    foreach_dimension()
        n.x.refine = n.x.prolongation = refine_injection;

    alpha.n = n;
    alpha.refine = alpha.prolongation = alpha_refine;
#endif
}


/*
The function below fills frag with the normal vector n, alpha, and the volume fraction of the
cell that is closest to the ghost cell. It also returns the coordinates of the fragment's midpoint 
and fills fluidCell with the cell center coordinates of the closest fluid cell.

It now also fills "bioff" with the indices of the closest interface w.r.t the ghost cell's stencil.

Note: we make a crude approximation that the closest interfacial point from the surrounding
cells to the ghost cell is which ever interfacial mid point/centroid is closest. In practice, it has worked
adequately, but this can be improved.

TODO: Clean up and streamline function.
*/

coord closest_interface (Point point, vector midPoints, scalar ibm, vector normals,
                         fragment * frag, coord * fluidCell, PointIBM * bioff)
{
    fragment temp_frag;
    coord temp_midPoint, temp_fluidCell = {0,0};
    coord n;
    PointIBM ptemp = {0,0,0};

    double min_distance = 1e6;

     for(int i = -1; i <= 1; i++) {
        double dx = midPoints.x[i] - x;
        double dy = midPoints.y[i] - y;
        double dz = midPoints.z[i] - z;
        if ((midPoints.x[i] || midPoints.y[i] || midPoints.z[i]) && 
             distance3D(dx, dy, dz) < min_distance) {
            temp_midPoint.x = midPoints.x[i];
            temp_midPoint.y = midPoints.y[i];
            temp_midPoint.z = midPoints.z[i];

            n.x = normals.x[i]; n.y = normals.y[i]; n.z = normals.z[i];

            fill_fragment (ibm[i], n, &temp_frag);
            temp_fluidCell.x = i*Delta + x;
            temp_fluidCell.y = y;
            temp_fluidCell.z = z;
            min_distance = distance3D(dx, dy, dz);
            ptemp = (PointIBM){i,0,0};
        }
     }

     for(int j = -1; j <= 1; j++) {
        double dx = midPoints.x[0,j] - x;
        double dy = midPoints.y[0,j] - y;
        double dz = midPoints.z[0,j] - z;
        if ((midPoints.x[0,j] || midPoints.y[0,j] || midPoints.z[0,j]) &&
             distance3D(dx, dy, dz) < min_distance) {
            temp_midPoint.x = midPoints.x[0,j];
            temp_midPoint.y = midPoints.y[0,j];
            temp_midPoint.z = midPoints.z[0,j];

            n.x = normals.x[0,j]; n.y = normals.y[0,j]; n.z = normals.z[0,j];

            fill_fragment (ibm[0,j], n, &temp_frag);
            temp_fluidCell.x = x;
            temp_fluidCell.y = j*Delta + y;
            temp_fluidCell.z = z;
            min_distance = distance3D(dx, dy, dz);
            ptemp = (PointIBM){0,j,0};
        }
     }
#if dimension == 3
     for(int k = -1; k <= 1; k++) {
        double dx = midPoints.x[0,0,k] - x;
        double dy = midPoints.y[0,0,k] - y;
        double dz = midPoints.z[0,0,k] - z;
        if ((midPoints.x[0,0,k] || midPoints.y[0,0,k] || midPoints.z[0,0,k]) &&
             distance3D(dx, dy, dz) < min_distance) {
            temp_midPoint.x = midPoints.x[0,0,k];
            temp_midPoint.y = midPoints.y[0,0,k];
            temp_midPoint.z = midPoints.z[0,0,k];

            n.x = normals.x[0,0,k]; n.y = normals.y[0,0,k]; n.z = normals.z[0,0,k];

            fill_fragment (ibm[0,0,k], n, &temp_frag);
            temp_fluidCell.x = x;
            temp_fluidCell.y = y;
            temp_fluidCell.z = k*Delta + z;
            min_distance = distance3D(dx, dy, dz);
            ptemp = (PointIBM){0,0,k};
        }
     }
#endif

    *fluidCell = temp_fluidCell;
    *frag = temp_frag;
    *bioff = ptemp;

    return temp_midPoint;
}


/*
The function below returns the boundary intercept coordinate given a fragment,
fluid cell coordinates, and volume fraction field (ibm).

TODO: Show derivation.
TODO: Handle degenerative case when boundary intercept is outside of cell.
*/

coord boundary_int (Point point, fragment frag, coord fluidCell, scalar ibm)
{
    double mag = distance3D(frag.n.x, frag.n.y, frag.n.z) + SEPS;
    coord n = frag.n, ghostCell = {x,y,z};

    normalize2(&n);

    double offset = 0;
    offset += n.x * -sign2(fluidCell.x - x);
    offset += n.y * -sign2(fluidCell.y - y);
    offset += n.z * -sign2(fluidCell.z - z);
    coord boundaryInt = {(-frag.alpha / mag - offset) * n.x, // is - correct? or should be positive alpha?
                         (-frag.alpha / mag - offset) * n.y,
                         (-frag.alpha / mag - offset) * n.z};

    foreach_dimension()
        boundaryInt.x = ghostCell.x + boundaryInt.x*Delta;

    return boundaryInt;
}

coord direction_vector (coord p0, coord p1)
{
    return (coord){p1.x - p0.x, p1.y - p0.y, p1.z - p0.z};
}


/*
image_point takes in the coordinates of the boundary intercept and ghost cell and 
returns the coordinates of the image point.
*/

coord image_point (coord boundaryInt, coord ghostCell)
{
     double dx = boundaryInt.x - ghostCell.x;
     double dy = boundaryInt.y - ghostCell.y;
     double dz = boundaryInt.z - ghostCell.z;

     coord imagePoint = {ghostCell.x + 2*dx, ghostCell.y + 2*dy, ghostCell.z + 2*dz};

     return imagePoint;
}


/*
Similarly to image_point, fresh_image_point calculates the imagepoint for "fresh cells"
*/

coord fresh_image_point (coord boundaryInt, coord freshCell)
{
     double dx = freshCell.x - boundaryInt.x;
     double dy = freshCell.y - boundaryInt.y;
     double dz = freshCell.z - boundaryInt.z;

     coord imagePoint = {freshCell.x + dx, freshCell.y + dy, freshCell.z + dz};

     return imagePoint;
}


/*
borders_boundary checks to see if a cell is adjacent to a domain boundary. 
Returns true if it has one or more neighbors that are inside the boundary.

Optionally, the user can provide two integer pointers (useri and userj) to be
filled with the indexs of the boundary neighbor.
*/

bool borders_boundary (Point point, int * useri = NULL, int * userj = NULL, int * userk = NULL)
{
#if TREE
    // Look at directly adjacent neighbors (4 in 2D)
    for (int d = 0; d < dimension; d++) {
	    for (int kk = -1; kk <= 1; kk += 2) {
            int _i = 0, _j = 0, _k = 0;
	        if (d == 0)
                _i = kk; 
            else if (d == 1)
                _j = kk; 
            else if (d == 2)
                _k = kk;
            // check to see if neighboring cell is inside boundary
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
/*
borders_boundary checks to see if a cell is adjacent to a MPI boundary. 
Returns true if it has one or more neighbors that are inside the boundary.
*/

bool borders_mpi_boundary (Point point)
{
#if TREE 
    // Look at directly adjacent neighbors (4 in 2D, 6 in 3D)
    for (int d = 0; d < dimension; d++) {
	    for (int kk = -1; kk <= 1; kk += 2) {
            int i = 0, j = 0, k = 0;
	        if (d == 0)
                i = kk; 
            else if (d == 1)
                j = kk; 
            else if (d == 2)
                k = kk;
            // check to see if neighboring cell is outside the current MPI rank/processor
            if (neighbor(-i,-j,-k).pid != pid())
                return true;
            (void) k;
        }
    }
#endif
    return false;
}
#endif


/*
This function extrapolates a scalar field p to a specified point given its coordinates,
normal vector n, and volume fraction field s.

TODO: cleaner 3D implementation.
*/

double extrapolate_scalar (Point point, scalar s, coord interpolatePoint, coord n, scalar p)
{
#if dimension == 2
    double weight[5][5] = {0};
    double weightSum = 0.;
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            if (s[i,j] > 0.5) {

                coord cellCenter = {x + Delta*i,y + Delta*j}, d;
                foreach_dimension()
                    d.x = interpolatePoint.x - cellCenter.x;

                double distanceMag = distance (d.x, d.y);
                double normalProjection = (n.x * d.x) + (n.y * d.y);

                weight[i][j] = sq(distanceMag) * fabs(normalProjection);

                weightSum += weight[i][j];
            }
            else
                weight[i][j] = 0.;
        }
    }

    double interpolatedScalar = 0;

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            interpolatedScalar += (weight[i][j]/(weightSum + SEPS)) * p[i,j];
        }
    }
#else
    double weight[5][5][5] = {0};
    double weightSum = 0.;
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            for (int k = -2; k <= 2; k++) {
                if (s[i,j,k] > 0.5) {

                    coord cellCenter = {x + Delta*i, y + Delta*j, z + Delta*k}, d;
                    foreach_dimension()
                        d.x = interpolatePoint.x - cellCenter.x;

                    double distanceMag = distance3D (d.x, d.y, d.z);
                    double normalProjection = (n.x * d.x) + (n.y * d.y) + (n.z * d.z);

                    weight[i+2][j+2][k+2] = sq(distanceMag) * fabs(normalProjection);

                    weightSum += weight[i+2][j+2][k+2];
                }
                else
                   weight[i+2][j+2][k+2] = 0.;
            }
        }
    }

    double interpolatedScalar = 0;

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            for (int k = -2; k <= 2; k++) {
                interpolatedScalar += (weight[i+2][j+2][k+2]/(weightSum + SEPS)) * p[i,j,k];
            }
        }
    }
#endif
    return interpolatedScalar;
}


/*
gauss_elim performs *in place* transformation to the provided augmented matrix 
(meaning it is changed w/o making a copy) and fills coeff with the solved linear system.

***Courtesy of ChatGPT***

TODO: Extend to handle higher-order interpolation schemes, i.e. larger matrices.
*/

void gauss_elim(int m, int n, double matrix[m][n], double sol[m])
{
    // Forward elimination
    for (int i = 0; i < m; i++) {

        // 1. Partial pivot: find row with largest pivot in column i
        int max_row = i;
        for (int r = i + 1; r < m; r++) {
            if (fabs(matrix[r][i]) > fabs(matrix[max_row][i])) {
                max_row = r;
            }
        }

        // 2. Swap current row i with max_row if needed
        if (max_row != i) {
            for (int c = 0; c < n; c++) {
                double temp = matrix[i][c];
                matrix[i][c]  = matrix[max_row][c];
                matrix[max_row][c] = temp;
            }
        }

        // 3. Make sure our pivot is non‐zero (or not too close to zero)
        if (fabs(matrix[i][i]) < 1e-12) {
            fprintf(stderr, "ERROR: Pivot is zero (matrix is singular or nearly singular)\n");
            return;
        }

        // 4. Eliminate all rows below row i
        for (int r = i + 1; r < m; r++) {
            double factor = matrix[r][i] / matrix[i][i];

            for (int c = i; c < n; c++) {
                matrix[r][c] -= factor * matrix[i][c];
            }
        }
    }

    // Back‐substitution
    for (int i = m - 1; i >= 0; i--) {
        // Start with the RHS of the augmented matrix
        sol[i] = matrix[i][n - 1];

        // Subtract the known terms from columns to the right
        for (int c = i + 1; c < m; c++) {
            sol[i] -= matrix[i][c] * sol[c];
        }

        // Divide by the diagonal element
        sol[i] /= matrix[i][i];
    }
}


/*
image_offsets fills integers xOffset and yOffset with the index of the cell containing
the given image point w.r.t the ghost cell's stencil.

e.g., if an image point is 2 cells to the left of it's respective ghost cell, then
xOffset = -2 and yOffset = 0.

*/

int image_offsets (Point point, coord imagePoint, int *xOffset, int *yOffset, int *zOffset = NULL)
{
    coord ghostCell = {x,y,z};
    int offset_x = 0, offset_y = 0, offset_z = 0;
    foreach_dimension() {
        double d = fabs(imagePoint.x - ghostCell.x) / Delta;
        int dsign = sign (imagePoint.x - ghostCell.x);
        if (d >= 1.5) {
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


/*
fluid_only checks to see if a point being used for interpolation is inside the solid
domain. If it is, the coordinates of that point is moved the a point on the interface,
which in this case is the interfacial midpoint. The velocity for this point is then
changed to the imposed boundary condition.

TODO: add check for completely full cells (ibm = 0)
TODO: Not able to do neumann for u right now (too many unknowns in interpolation)?
TODO: allow for navier-slip condition!
*/

extern vector u;

void fluid_only (Point point, int xx, int yy, int zz, int i, int j, int k, 
                 coord * pTemp, coord * velocity, vector midPoints,
                 int bOffset_X, int bOffset_Y, int bOffset_Z, 
                 vector normals, coord imagePoint)
{
    int off_x = xx + i, off_y = yy + j, off_z = zz + k;
    if (ibm[off_x,off_y,off_z] <= 0.5 && ibm[off_x,off_y,off_z] > 0.) {
        pTemp->x = midPoints.x[off_x,off_y,off_z];
        pTemp->y = midPoints.y[off_x,off_y,off_z];
        pTemp->z = midPoints.z[off_x,off_y,off_z];

        if (local_bc_coordinates) {
            coord gcvelocity = {u.x[off_x, off_y, off_z],
                                u.y[off_x, off_y, off_z],
                                u.z[off_x, off_y, off_z]};
            coord n = {normals.x[off_x, off_y, off_z],
                       normals.y[off_x, off_y, off_z],
                       normals.z[off_x, off_y, off_z]};
            coord t1, t2;
            normal_and_tangents (&n, &t1, &t2);
            coord gcprojVelocity = {dot_product(gcvelocity, n),
                                    dot_product(gcvelocity, t1),
                                    dot_product(gcvelocity, t2)};

            coord projVelocity = {0,0,0};
            //coord ghostCell = {x + off_x*Delta, y + off_y*Delta, z + off_z*Delta};
            coord boundp = {midPoints.x[off_x,off_y,off_z], 
                            midPoints.y[off_x,off_y,off_z], 
                            midPoints.z[off_x,off_y,off_z]};

            coord d = {-n.x, -n.y, -n.z}; // inward pointing normal on solid interface

            foreach_dimension() {
                if (bOffset_X == off_x) 
                    pTemp->x += bOffset_X * Delta;
                
                bool bctype[2] = {false, false};
                // TODO: this uses the base ghost cells B.C, when it should use
                //       the one that the interpolation point is inside of! use a foreach_point?
                double vb = u.x.boundary[immersed] (point, point, u.x, bctype);
                bool dirichlet = bctype[0], nslip = bctype[1];
                if (dirichlet) {
                    if (nslip) {
                        projVelocity.x = 0; // navier bc doesn't work with moving bodies for now
                        foreach_dimension() // change point
                            pTemp->x = boundp.x  + vb*d.x; // vb = slip length
                    }
                    else {
                        projVelocity.x = vb;
                    }
                }
                else { // neumann
                    projVelocity.x = gcprojVelocity.x;
                    //foreach_dimension()
                    //    pTemp->x = ghostCell.x;
                }
            }

            double gcn = projVelocity.x, gct1 = projVelocity.y, gct2 = projVelocity.z;
            foreach_dimension()
                velocity->x = gcn*n.x + gct1*t1.x + gct2*t2.x;
        }
        else { // !local_bc_coordinates, i.e. use n ≡ x, t ≡ y, and r ≡ z.
            foreach_dimension() {
                bool dirichlet = true;
                double vb = u.x.boundary[immersed] (point, point, u.x, &dirichlet);
                if (dirichlet)
                    velocity->x = vb;
                else
                    velocity->x = u.x[xx + i, yy + j, zz + k];
                }
        }
    }
    (void) off_z; // to prevent unused variable warning
}


/*
The function below uses interpolation to find the velocity at the image point and
returns it given a vector field, u, the coordinates of the image point, and a field
containing all interfacial midpoints.

TODO: Streamline and clean-up code?
TODO: Extend to handle higher-order interpolation schemes, i.e. larger matrices.
*/

coord image_velocity (Point point, vector u, coord imagePoint, vector midPoints, vector normals)
{
    
    int boundaryOffsetX = 0, boundaryOffsetY = 0, boundaryOffsetZ = 0;
    borders_boundary (point, &boundaryOffsetX, &boundaryOffsetY, &boundaryOffsetZ);
    
    int xOffset = 0, yOffset = 0, zOffset = 0;
    image_offsets (point, imagePoint, &xOffset, &yOffset, &zOffset);
    
    assert (abs(xOffset) <= 2 && abs(yOffset) <= 2 && abs(zOffset) <= 2);

    coord imageCell = {x + Delta*xOffset, y + Delta*yOffset, z + Delta*zOffset};
    
    int i = sign(imagePoint.x - imageCell.x);
    int j = sign(imagePoint.y - imageCell.y);
    int k = sign(imagePoint.z - imageCell.z);

    int xx = xOffset, yy = yOffset, zz = zOffset;

    coord velocity[(int)pow(2, dimension)]; // 4 in 2D, 8 in 3D
    velocity[0].x = u.x[xx,yy,zz];
    velocity[1].x = u.x[xx+i,yy,zz];
    velocity[2].x = u.x[xx+i,yy+j,zz];
    velocity[3].x = u.x[xx,yy+j,zz];

    velocity[0].y = u.y[xx,yy,zz];
    velocity[1].y = u.y[xx+i,yy,zz];
    velocity[2].y = u.y[xx+i,yy+j,zz];
    velocity[3].y = u.y[xx,yy+j,zz];

#if dimension == 3
    velocity[4].x = u.x[xx,yy,zz+k];
    velocity[5].x = u.x[xx+i,yy,zz+k];
    velocity[6].x = u.x[xx+i,yy+j,zz+k];
    velocity[7].x = u.x[xx,yy+j,zz+k];

    velocity[4].y = u.y[xx,yy,zz+k];
    velocity[5].y = u.y[xx+i,yy,zz+k];
    velocity[6].y = u.y[xx+i,yy+j,zz+k];
    velocity[7].y = u.y[xx,yy+j,zz+k];

    velocity[0].z = u.z[xx,yy,zz];
    velocity[1].z = u.z[xx+i,yy,zz];
    velocity[2].z = u.z[xx+i,yy+j,zz];
    velocity[3].z = u.z[xx,yy+j,zz];
    velocity[4].z = u.z[xx,yy,zz+k];
    velocity[5].z = u.z[xx+i,yy,zz+k];
    velocity[6].z = u.z[xx+i,yy+j,zz+k];
    velocity[7].z = u.z[xx,yy+j,zz+k];
#endif

    coord p0 = {imageCell.x, imageCell.y, imageCell.z};
    coord p1 = {imageCell.x + i*Delta, imageCell.y, imageCell.z};
    coord p2 = {imageCell.x + i*Delta, imageCell.y + j*Delta, imageCell.z};
    coord p3 = {imageCell.x, imageCell.y + j*Delta, imageCell.z};
#if dimension == 3
    coord p4 = {imageCell.x, imageCell.y, imageCell.z + k*Delta};
    coord p5 = {imageCell.x + i*Delta, imageCell.y, imageCell.z + k*Delta};
    coord p6 = {imageCell.x + i*Delta, imageCell.y + j*Delta, imageCell.z + k*Delta};
    coord p7 = {imageCell.x, imageCell.y + j*Delta, imageCell.z + k*Delta};
#endif
   
    // make sure all points are inside the fluid domain ...
    // if not, change their coordinates to a point on the interface
    fluid_only (point, xx, yy, zz, 0, 0, 0, &p0, &velocity[0], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);

    fluid_only (point, xx, yy, zz, i, 0, 0, &p1, &velocity[1], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);

    fluid_only (point, xx, yy, zz, i, j, 0, &p2, &velocity[2], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);

    fluid_only (point, xx, yy, zz, 0, j, 0, &p3, &velocity[3], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);
#if dimension == 3
    fluid_only (point, xx, yy, zz, 0, 0, k, &p4, &velocity[4], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);

    fluid_only (point, xx, yy, zz, i, 0, k, &p5, &velocity[5], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);

    fluid_only (point, xx, yy, zz, i, j, k, &p6, &velocity[6], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);

    fluid_only (point, xx, yy, zz, 0, j, k, &p7, &velocity[7], midPoints, 
                boundaryOffsetX, boundaryOffsetY, boundaryOffsetZ, normals, imagePoint);
#endif

#if dimension == 2
    double vanderVelo_x[4][5] = {
        {p0.x*p0.y, p0.x, p0.y, 1, velocity[0].x},
        {p1.x*p1.y, p1.x, p1.y, 1, velocity[1].x},
        {p2.x*p2.y, p2.x, p2.y, 1, velocity[2].x},
        {p3.x*p3.y, p3.x, p3.y, 1, velocity[3].x},
    };

    double vanderVelo_y[4][5] = {
        {p0.x*p0.y, p0.x, p0.y, 1, velocity[0].y},
        {p1.x*p1.y, p1.x, p1.y, 1, velocity[1].y},
        {p2.x*p2.y, p2.x, p2.y, 1, velocity[2].y},
        {p3.x*p3.y, p3.x, p3.y, 1, velocity[3].y},
    };

    int m = 4, n = 5;
    double coeff_x[4], coeff_y[4];
#else // dimension = 3
    double vanderVelo_x[8][9] = {
        {p0.x*p0.y*p0.z, p0.x*p0.y, p0.x*p0.z, p0.y*p0.z, p0.x, p0.y, p0.z, 1, velocity[0].x},
        {p1.x*p1.y*p1.z, p1.x*p1.y, p1.x*p1.z, p1.y*p1.z, p1.x, p1.y, p1.z, 1, velocity[1].x},
        {p2.x*p2.y*p2.z, p2.x*p2.y, p2.x*p2.z, p2.y*p2.z, p2.x, p2.y, p2.z, 1, velocity[2].x},
        {p3.x*p3.y*p3.z, p3.x*p3.y, p3.x*p3.z, p3.y*p3.z, p3.x, p3.y, p3.z, 1, velocity[3].x},
        {p4.x*p4.y*p4.z, p4.x*p4.y, p4.x*p4.z, p4.y*p4.z, p4.x, p4.y, p4.z, 1, velocity[4].x},
        {p5.x*p5.y*p5.z, p5.x*p5.y, p5.x*p5.z, p5.y*p5.z, p5.x, p5.y, p5.z, 1, velocity[5].x},
        {p6.x*p6.y*p6.z, p6.x*p6.y, p6.x*p6.z, p6.y*p6.z, p6.x, p6.y, p6.z, 1, velocity[6].x},
        {p7.x*p7.y*p7.z, p7.x*p7.y, p7.x*p7.z, p7.y*p7.z, p7.x, p7.y, p7.z, 1, velocity[7].x},
    };

    double vanderVelo_y[8][9] = {
        {p0.x*p0.y*p0.z, p0.x*p0.y, p0.x*p0.z, p0.y*p0.z, p0.x, p0.y, p0.z, 1, velocity[0].y},
        {p1.x*p1.y*p1.z, p1.x*p1.y, p1.x*p1.z, p1.y*p1.z, p1.x, p1.y, p1.z, 1, velocity[1].y},
        {p2.x*p2.y*p2.z, p2.x*p2.y, p2.x*p2.z, p2.y*p2.z, p2.x, p2.y, p2.z, 1, velocity[2].y},
        {p3.x*p3.y*p3.z, p3.x*p3.y, p3.x*p3.z, p3.y*p3.z, p3.x, p3.y, p3.z, 1, velocity[3].y},
        {p4.x*p4.y*p4.z, p4.x*p4.y, p4.x*p4.z, p4.y*p4.z, p4.x, p4.y, p4.z, 1, velocity[4].y},
        {p5.x*p5.y*p5.z, p5.x*p5.y, p5.x*p5.z, p5.y*p5.z, p5.x, p5.y, p5.z, 1, velocity[5].y},
        {p6.x*p6.y*p6.z, p6.x*p6.y, p6.x*p6.z, p6.y*p6.z, p6.x, p6.y, p6.z, 1, velocity[6].y},
        {p7.x*p7.y*p7.z, p7.x*p7.y, p7.x*p7.z, p7.y*p7.z, p7.x, p7.y, p7.z, 1, velocity[7].y},
    };

    double vanderVelo_z[8][9] = {
        {p0.x*p0.y*p0.z, p0.x*p0.y, p0.x*p0.z, p0.y*p0.z, p0.x, p0.y, p0.z, 1, velocity[0].z},
        {p1.x*p1.y*p1.z, p1.x*p1.y, p1.x*p1.z, p1.y*p1.z, p1.x, p1.y, p1.z, 1, velocity[1].z},
        {p2.x*p2.y*p2.z, p2.x*p2.y, p2.x*p2.z, p2.y*p2.z, p2.x, p2.y, p2.z, 1, velocity[2].z},
        {p3.x*p3.y*p3.z, p3.x*p3.y, p3.x*p3.z, p3.y*p3.z, p3.x, p3.y, p3.z, 1, velocity[3].z},
        {p4.x*p4.y*p4.z, p4.x*p4.y, p4.x*p4.z, p4.y*p4.z, p4.x, p4.y, p4.z, 1, velocity[4].z},
        {p5.x*p5.y*p5.z, p5.x*p5.y, p5.x*p5.z, p5.y*p5.z, p5.x, p5.y, p5.z, 1, velocity[5].z},
        {p6.x*p6.y*p6.z, p6.x*p6.y, p6.x*p6.z, p6.y*p6.z, p6.x, p6.y, p6.z, 1, velocity[6].z},
        {p7.x*p7.y*p7.z, p7.x*p7.y, p7.x*p7.z, p7.y*p7.z, p7.x, p7.y, p7.z, 1, velocity[7].z},
    };

    int m = 8, n = 9;
    double coeff_x[8], coeff_y[8], coeff_z[8];
#endif

    foreach_dimension()
        gauss_elim (m, n, vanderVelo_x, coeff_x);

    coord temp_velo = {0,0,0};

#if dimension == 2
    temp_velo.x = coeff_x[0] * imagePoint.x * imagePoint.y +
                  coeff_x[1] * imagePoint.x +
                  coeff_x[2] * imagePoint.y +
                  coeff_x[3];

    temp_velo.y = coeff_y[0] * imagePoint.x * imagePoint.y +
                  coeff_y[1] * imagePoint.x +
                  coeff_y[2] * imagePoint.y +
                  coeff_y[3];
#else
    temp_velo.x = coeff_x[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                  coeff_x[1] * imagePoint.x * imagePoint.y +
                  coeff_x[2] * imagePoint.x * imagePoint.z +
                  coeff_x[3] * imagePoint.y * imagePoint.z +
                  coeff_x[4] * imagePoint.x +
                  coeff_x[5] * imagePoint.y +
                  coeff_x[6] * imagePoint.z +
                  coeff_x[7];

    temp_velo.y = coeff_y[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                  coeff_y[1] * imagePoint.x * imagePoint.y +
                  coeff_y[2] * imagePoint.x * imagePoint.z +
                  coeff_y[3] * imagePoint.y * imagePoint.z +
                  coeff_y[4] * imagePoint.x +
                  coeff_y[5] * imagePoint.y +
                  coeff_y[6] * imagePoint.z +
                  coeff_y[7];

    temp_velo.z = coeff_z[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                  coeff_z[1] * imagePoint.x * imagePoint.y +
                  coeff_z[2] * imagePoint.x * imagePoint.z +
                  coeff_z[3] * imagePoint.y * imagePoint.z +
                  coeff_z[4] * imagePoint.x +
                  coeff_z[5] * imagePoint.y +
                  coeff_z[6] * imagePoint.z +
                  coeff_z[7];
#endif
    
    (void) zz; (void) k;

    return temp_velo;

}


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

#define rows ((int)pow(2,dimension))
#define cols ((int)pow(2,dimension) + 1)

// takes in a row
// TODO: does the projected velocity hold true when using another ghost cell's interface coordinate system?
// e.g., when the left cell is a ghost cell, do we project u according to that cells n,t1,and t2? or keep it
// with the "home/center" ghost cell.

void fluid_only2 (Point point, const int n, double rmatrix[n],
                  PointIBM poff, PointIBM pnode, PointIBM pbound, 
                  char dir, coord * pcell, coord velocity, coord ipoint,
                  vector midPoints, vector normals, scalar alphas)
{
    // TODO: this uses the base ghost cells B.C, when it should use
    //       the one that the interpolation point is inside of! use a foreach_point?

    double val = 0;
    if      (dir == 'n') val = velocity.x;
    else if (dir == 't') val = velocity.y;
    else if (dir == 'r') val = velocity.z;

    int xx = poff.i + pnode.i, yy = poff.j + pnode.j, zz = poff.k + pnode.k;
    int type = 0;

    // a. Check to see if point is in a ghost cell, if so, move it to the interface
    //    and recalculate the node's value considering the immersed boundary condition.
    if (ibm[xx,yy,zz] <= 0.5 && ibm[xx,yy,zz] > 0.) {
        *pcell = (coord){midPoints.x[xx,yy,zz], midPoints.y[xx,yy,zz], midPoints.z[xx,yy,zz]};

        // move point more if cell is inside domain boundary
        if (xx == pbound.i) pcell->x += pbound.i*Delta;
        if (yy == pbound.j) pcell->y += pbound.j*Delta;
        if (zz == pbound.k) pcell->z += pbound.k*Delta;

        // b. check the bc type and get its value
        bool bctype[2] = {false, false};
        double bc = 0;

        if      (dir == 'n') bc = u.x.boundary[immersed] (point, point, u.x, bctype);
        else if (dir == 't') bc = u.y.boundary[immersed] (point, point, u.y, bctype);
#if dimension == 3
        else if (dir == 'r') bc = u.z.boundary[immersed] (point, point, u.z, bctype);
#endif
        bool dirichlet = bctype[0], navierslip = bctype[1];

        if (dirichlet && !navierslip) {
            val = bc;
            type = 0;
        }
        else if (dirichlet && navierslip) {
            assert(dir != 'n'); // navier-slip cannot be applied in the normal direction
            val = bc; // assumes stationary wall
            type = 2;
        }
        else { // neumann
            val = bc;
            type = 1;
        }
    }

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
        coord n = {normals.x[xx,yy,zz], normals.y[xx,yy,zz], normals.z[xx,yy,zz]};
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
        coord n = {normals.x[xx,yy,zz], normals.y[xx,yy,zz], normals.z[xx,yy,zz]};
        double mag = distance3D(n.x, n.y, n.z) + SEPS;
        normalize(&n);
        double alpha = alphas[xx,yy,zz];

        coord gc = {x + xx*Delta, y + yy*Delta, z + zz*Delta}, bi; // ghost cell, boundary intercept
        foreach_dimension()
            bi.x = alpha/mag * n.x;

        foreach_dimension()
            pcell->x = gc.x + bi.x*Delta;

        coord d = direction_vector(gc, *pcell);
        double dmag = distance3D(d.x, d.y, d.z); // distance from ghost cell to boundary intercept
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

    (void) zz; // to prevent unused variable warning
}


void get_interpolation_matrix (Point point, int m, int n, double matrix[m][n], char dir,
                               coord velo[m], PointIBM poff, PointIBM pnode, PointIBM pbound, 
                               coord ipoint, vector midPoints, vector normals, scalar alphas)
{
    // 4.a Calculate (global) coordinates and relative indices of interpolation cells
    coord pints[rows];
    PointIBM pnodes[rows];
    get_interpolation_points(point, m, pints, pnodes, poff, pnode);

    // 4.b Assemble matrix row by row
    for (int row = 0; row < m; ++row) {

        // If a cell for interpolating is a ghost cell, move the point to the
        // interface and change the velocity to the correct boundary condition
        fluid_only2(point, n, matrix[row], poff, pnodes[row], pbound, dir, 
                    &pints[row], velo[row], ipoint, midPoints, normals, alphas);
    }
}

coord image_velocity2 (Point point, vector u, coord imagePoint, PointIBM bioff, 
                       vector midPoints, vector normals, scalar alphas)
{
    // 1. Calculate offsets 
    int boffx = 0, boffy = 0, boffz = 0; // boundary offsets
    borders_boundary (point, &boffx, &boffy, &boffz);
    
    int xOffset = 0, yOffset = 0, zOffset = 0;
    image_offsets (point, imagePoint, &xOffset, &yOffset, &zOffset);
    
    assert (abs(xOffset) <= 2 && abs(yOffset) <= 2 && abs(zOffset) <= 2);

    coord imageCell = {x + Delta*xOffset, y + Delta*yOffset, z + Delta*zOffset};
    
    int i = sign(imagePoint.x - imageCell.x);
    int j = sign(imagePoint.y - imageCell.y);
    int k = sign(imagePoint.z - imageCell.z);

    int xx = xOffset, yy = yOffset, zz = zOffset;

    // 2. Grab velocity from cells used for interpolation 
    //TODO: condense this and maybe move to separate function
    coord velocity[rows];
    velocity[0].x = u.x[xx,yy,zz];
    velocity[1].x = u.x[xx+i,yy,zz];
    velocity[2].x = u.x[xx+i,yy+j,zz];
    velocity[3].x = u.x[xx,yy+j,zz];

    velocity[0].y = u.y[xx,yy,zz];
    velocity[1].y = u.y[xx+i,yy,zz];
    velocity[2].y = u.y[xx+i,yy+j,zz];
    velocity[3].y = u.y[xx,yy+j,zz];

#if dimension == 3
    velocity[4].x = u.x[xx,yy,zz+k];
    velocity[5].x = u.x[xx+i,yy,zz+k];
    velocity[6].x = u.x[xx+i,yy+j,zz+k];
    velocity[7].x = u.x[xx,yy+j,zz+k];

    velocity[4].y = u.y[xx,yy,zz+k];
    velocity[5].y = u.y[xx+i,yy,zz+k];
    velocity[6].y = u.y[xx+i,yy+j,zz+k];
    velocity[7].y = u.y[xx,yy+j,zz+k];

    velocity[0].z = u.z[xx,yy,zz];
    velocity[1].z = u.z[xx+i,yy,zz];
    velocity[2].z = u.z[xx+i,yy+j,zz];
    velocity[3].z = u.z[xx,yy+j,zz];
    velocity[4].z = u.z[xx,yy,zz+k];
    velocity[5].z = u.z[xx+i,yy,zz+k];
    velocity[6].z = u.z[xx+i,yy+j,zz+k];
    velocity[7].z = u.z[xx,yy+j,zz+k];
#endif

    // 3. Project velocity to normal and tangent(s) direction
    coord n = {normals.x[bioff.i,bioff.j,bioff.k], 
               normals.y[bioff.i,bioff.j,bioff.k], 
               normals.z[bioff.i,bioff.j,bioff.k]}, t1, t2;
    normal_and_tangents (&n, &t1, &t2);

    for (int i = 0; i < rows; ++i) {
        coord projvelo = {dot_product(velocity[i], n),
                          dot_product(velocity[i], t1),
                          dot_product(velocity[i], t2)};
        velocity[i] = (coord){projvelo.x, projvelo.y, projvelo.z};
    }

    // 4. Assemble interpolation matrices
    double veloMatrix_x[rows][cols];
    get_interpolation_matrix(point, rows, cols, veloMatrix_x, 'n', velocity, 
                            (PointIBM){xx,yy,zz}, (PointIBM){i,j,k}, 
                            (PointIBM){boffx,boffy,boffz}, imagePoint, 
                             midPoints, normals, alphas);
    double veloMatrix_y[rows][cols];
    get_interpolation_matrix(point, rows, cols, veloMatrix_y, 't', velocity,
                            (PointIBM){xx,yy,zz}, (PointIBM){i,j,k}, 
                            (PointIBM){boffx,boffy,boffz}, imagePoint, 
                             midPoints, normals, alphas);

    double coeff_x[rows], coeff_y[rows];

#if dimension == 3
    double veloMatrix_z[rows][cols];
    get_interpolation_matrix(point, rows, cols, veloMatrix_z, 'r', velocity,
                            (PointIBM){xx,yy,zz}, (PointIBM){i,j,k}, 
                            (PointIBM){boffx,boffy,boffz}, imagePoint, 
                             midPoints, normals, alphas);
    double coeff_z[rows];
#endif

    // 5. Solve the system of linear equations to get the interpolating coefficients
    foreach_dimension()
        gauss_elim (rows, cols, veloMatrix_x, coeff_x);

    // 6. Calculate the image velocity
    coord imageVelo = {0,0,0};
#if dimension == 2
    imageVelo.x = coeff_x[0] * imagePoint.x * imagePoint.y +
                  coeff_x[1] * imagePoint.x +
                  coeff_x[2] * imagePoint.y +
                  coeff_x[3];
 
    imageVelo.y = coeff_y[0] * imagePoint.x * imagePoint.y +
                  coeff_y[1] * imagePoint.x +
                  coeff_y[2] * imagePoint.y +
                  coeff_y[3];
#else
    imageVelo.x = coeff_x[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                  coeff_x[1] * imagePoint.x * imagePoint.y +
                  coeff_x[2] * imagePoint.x * imagePoint.z +
                  coeff_x[3] * imagePoint.y * imagePoint.z +
                  coeff_x[4] * imagePoint.x +
                  coeff_x[5] * imagePoint.y +
                  coeff_x[6] * imagePoint.z +
                  coeff_x[7];

    imageVelo.y = coeff_y[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                  coeff_y[1] * imagePoint.x * imagePoint.y +
                  coeff_y[2] * imagePoint.x * imagePoint.z +
                  coeff_y[3] * imagePoint.y * imagePoint.z +
                  coeff_y[4] * imagePoint.x +
                  coeff_y[5] * imagePoint.y +
                  coeff_y[6] * imagePoint.z +
                  coeff_y[7];

    imageVelo.z = coeff_z[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                  coeff_z[1] * imagePoint.x * imagePoint.y +
                  coeff_z[2] * imagePoint.x * imagePoint.z +
                  coeff_z[3] * imagePoint.y * imagePoint.z +
                  coeff_z[4] * imagePoint.x +
                  coeff_z[5] * imagePoint.y +
                  coeff_z[6] * imagePoint.z +
                  coeff_z[7];
#endif

    // 7. Project the velocity back to the cartesian coordinate system
    // TODO: is this necessary if we project it back when calculating the gc value?
    double iux = imageVelo.x, iuy = imageVelo.y, iuz = imageVelo.z;
    foreach_dimension() 
        imageVelo.x = iux*n.x + iuy*t1.x + iuz*t2.x;
    
    (void) zz; (void) k;

    return imageVelo;
}

#undef rows
#undef cols

/*
The function below uses interpolation to find the velocity at the image point and
returns it given a vector field, u, the coordinates of the image point, and an array
containing all boundary intercept points.

TODO: Streamline and clean-up code.
TODO: Extend to handle higher-order interpolation schemes, i.e. larger matrices.
TODO: Combine with image_velocity to have only 1 function. Or generalize functions
      to handle any vector or scalar?
TODO: 
*/

double image_pressure (Point point, scalar p, coord imagePoint) 
{

    int xOffset = 0, yOffset = 0, zOffset;
    image_offsets (point, imagePoint, &xOffset, &yOffset, &zOffset);
    
    assert (abs(xOffset) <= 2 && abs(yOffset) <= 2 && abs(zOffset) <= 2);

    coord imageCell = {x + Delta * xOffset, y + Delta * yOffset, z + Delta * zOffset};

    int i = sign(imagePoint.x - imageCell.x);
    int j = sign(imagePoint.y - imageCell.y);
    int k = sign(imagePoint.z - imageCell.z);

    int xx = xOffset, yy = yOffset, zz = zOffset;

    double pressure[(int)pow(2, dimension)]; // 4 in 2D, 8 in 3D
    pressure[0] = p[xx,yy,zz];
    pressure[1] = p[xx+i,yy,zz];
    pressure[2] = p[xx+i,yy+j,zz];
    pressure[3] = p[xx,yy+j,zz];
#if dimension == 3
    pressure[4] = p[xx,yy,zz+k];
    pressure[5] = p[xx+i,yy,zz+k];
    pressure[6] = p[xx+i,yy+j,zz+k];
    pressure[7] = p[xx,yy+j,zz+k];
#endif

    coord p0 = {imageCell.x, imageCell.y, imageCell.z};
    coord p1 = {imageCell.x + i*Delta, imageCell.y, imageCell.z};
    coord p2 = {imageCell.x + i*Delta, imageCell.y + j*Delta, imageCell.z};
    coord p3 = {imageCell.x, imageCell.y + j*Delta, imageCell.z};
#if dimension == 3
    coord p4 = {imageCell.x, imageCell.y, imageCell.z + k*Delta};
    coord p5 = {imageCell.x + i*Delta, imageCell.y, imageCell.z + k*Delta};
    coord p6 = {imageCell.x + i*Delta, imageCell.y + j*Delta, imageCell.z + k*Delta};
    coord p7 = {imageCell.x, imageCell.y + j*Delta, imageCell.z + k*Delta};
#endif

#if dimension == 2
    double vanderPressure[4][5] = {
        {p0.x*p0.y, p0.x, p0.y, 1, pressure[0]},
        {p1.x*p1.y, p1.x, p1.y, 1, pressure[1]},
        {p2.x*p2.y, p2.x, p2.y, 1, pressure[2]},
        {p3.x*p3.y, p3.x, p3.y, 1, pressure[3]},
    };

    int m = 4, n = 5;
    double coeff[4];
#else
    double vanderPressure[8][9] = {
        {p0.x*p0.y*p0.z, p0.x*p0.y, p0.x*p0.z, p0.y*p0.z, p0.x, p0.y, p0.z, 1, pressure[0]},
        {p1.x*p1.y*p1.z, p1.x*p1.y, p1.x*p1.z, p1.y*p1.z, p1.x, p1.y, p1.z, 1, pressure[1]},
        {p2.x*p2.y*p2.z, p2.x*p2.y, p2.x*p2.z, p2.y*p2.z, p2.x, p2.y, p2.z, 1, pressure[2]},
        {p3.x*p3.y*p3.z, p3.x*p3.y, p3.x*p3.z, p3.y*p3.z, p3.x, p3.y, p3.z, 1, pressure[3]},
        {p4.x*p4.y*p4.z, p4.x*p4.y, p4.x*p4.z, p4.y*p4.z, p4.x, p4.y, p4.z, 1, pressure[4]},
        {p5.x*p5.y*p5.z, p5.x*p5.y, p5.x*p5.z, p5.y*p5.z, p5.x, p5.y, p5.z, 1, pressure[5]},
        {p6.x*p6.y*p6.z, p6.x*p6.y, p6.x*p6.z, p6.y*p6.z, p6.x, p6.y, p6.z, 1, pressure[6]},
        {p7.x*p7.y*p7.z, p7.x*p7.y, p7.x*p7.z, p7.y*p7.z, p7.x, p7.y, p7.z, 1, pressure[7]},
    };

    int m = 8, n = 9;
    double coeff[8];
#endif

    gauss_elim (m, n, vanderPressure, coeff);

#if dimension == 2
    double temp_pressure = coeff[0] * imagePoint.x*imagePoint.y +
                           coeff[1] * imagePoint.x +
                           coeff[2] * imagePoint.y +
                           coeff[3];
#else
    double temp_pressure = coeff[0] * imagePoint.x * imagePoint.y * imagePoint.z +
                           coeff[1] * imagePoint.x * imagePoint.y +
                           coeff[2] * imagePoint.x * imagePoint.z +
                           coeff[3] * imagePoint.y * imagePoint.z +
                           coeff[4] * imagePoint.x +
                           coeff[5] * imagePoint.y +
                           coeff[6] * imagePoint.z +
                           coeff[7];
#endif

    (void) zz; (void) k;

    return temp_pressure;
}


/*
This macro is for refining a scalar field s and its face vector field sf to make
sure its interface is at the specified max level. The user specifies a level set 
function in *expr* and also the maximum number of iterations.

max_i = 100 to 1000 seems to be enough.

Note, the user should call the solid function again after using this macro to
initialize it on the newly refined field.
*/

#define initial_refine(s, sf, expr, maxLevel, minLevel, max_i)    \
    astats ss;                                                    \
    int ic = 0;                                                   \
    do {                                                          \
        ic++;                                                     \
        solid (s, sf, expr);                                      \
        ss = adapt_wavelet ({s, sf}, (double[]) {1.e-30, 1.e-30}, \
                maxlevel = maxLevel, minlevel = minLevel);        \
    } while ((ss.nf || ss.nc) && ic < max_i);                     \
    refine (s[] > 0 && s[] < 1 && level < maxLevel);



/*
Again, this function is taken from embed.h. It removes any cell with inconsistent 
volume/surface fractions.
*/

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
   
    foreach_face()
      if (s.x[] && ((!c[] || !c[-1]) || s.x[] < smin))
	s.x[] = 0.;

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


/*
ibm_geometry is taken from embed's embed_geometry. Given a point containing a
interface fragment, the function returns the area of the fragment and fills p and
n with the coordinates of the interfacial mid point and normal vector, respectively.

Optionally, the user can provide a pointer 'alpha' to be filled with the fragments
alpha value.

Note, n is normalized differently from the MYC's way. Using the interface_normal
function from fractions.h, |n.x| + |n.y| = 1. Here, n is normalized using the vectors
magnitude, i.e. sqrt( sq(n.x) + sq(n.y) ) = 1.

Note area must be multiplied by the cell length (or area in 3D) to get in terms of
phyiscal "units". Also, the midpoint, p, is in the cell's local coordinate system.
*/

static inline
//double ibm_geometry (Point point, coord * p, coord * n, double * alphau = NULL)
double ibm_geometry (Point point, coord * p, coord * n)
{
    *n = facet_normal (point, ibm, ibmf);
    double alpha = plane_alpha (ibm[], *n);

    //if (alphau != NULL)
    //    *alphau = alpha;

    double area = plane_area_center (*n, alpha, p);
    foreach_dimension()
        n->x *= -1;
    normalize (n);

    return area;
}


static inline
double ibm0_geometry (Point point, coord * p, coord * n, scalar ibm1, face vector ibmf1)
{
    *n = facet_normal (point, ibm1, ibmf1);
    double alpha = plane_alpha (ibm1[], *n);
    double area = plane_area_center (*n, alpha, p);
    foreach_dimension()
        n->x *= -1;
    normalize (n);

    return area;
}


void normalize_norm (coord n, coord * newn)
{
    double norm = fabs(n.y) + fabs(n.x) + fabs(n.z);
    coord nNorm = {n.x/norm, n.y/norm, n.z/norm};
    *newn = nNorm;
}


/*
borders_ghost_x checks to see if a given cell shares an x face with a ghost
cell. If it does, the function returns the ghost cell's index w.r.t the given stencil.

The function is automatically changed to handle the y direction (top and bottom faces)
via the foreach_dimension() operator.
*/

foreach_dimension()
int borders_ghost_x (Point point, scalar ibm)
{
    for (int i = -1; i <= 1; i += 2) {
        if (ibm[i] < 0.5 && ibm[i] > 0) {
            return i;
        }
    }
    return 0;
}


/*
Typically, the default functions that calculates n reguires that the "point" 
be the cell in which we determine n. This function allows us to calculate
a neighbor's normal vector by accessing it's 5x5 stencil.

Note, since this uses face fractions, we can only calculate n for neighboring cells
within its 3x3 stencil.

TODO: verify normalization order is correct (use normalize() after traditional normalzation?)
TODO: 3D implementation
*/

coord offset_normal (Point point, face vector sf, int xoffset, int yoffset)
{
    assert (abs(xoffset) <= 1 && abs(yoffset) <= 1); // to avoid out-of-bounds access
                                                     // from 5x5 stencil
    double nx = sf.x[xoffset,yoffset] - sf.x[xoffset + 1,yoffset];
    double ny = sf.y[xoffset,yoffset] - sf.y[xoffset,yoffset + 1];

    double mag = distance (nx, ny);

    coord n;
    if (mag == 0) {
        foreach_dimension() {
            n.x = 1./dimension;
        }
    }
    else {
        n.x = nx / mag; n.y = ny / mag;
    }
    return n;
}


/*
ghost_fluxes returns the fluxes of a ghost cell (small cut-cell) that is to
be redistributed or virtually merged with non-ghost cell neighbors.

Note: THIS IS CURRENTLY NOT USED

TODO: clean up and streamline code
TODO: 3D implementation
TODO: what if ghost cell has two fluid cell neighbors in one direction
      (left and right or top and bottom)?
*/

#if 0
coord ghost_fluxes (Point point, scalar ibm, face vector ibmf, face vector uf)
{
    int xindex = is_mostly_solid (ibm, 0)? 0: borders_ghost_x (point, ibm);
    int yindex = is_mostly_solid (ibm, 0)? 0: borders_ghost_y (point, ibm);
#if 0
    fprintf (stderr, "### New Cell ###\n");
    fprintf (stderr, "|| xindex=%d yindex=%d\n", xindex, yindex);
#endif
    assert (abs(xindex) <= 1 && abs(yindex) <= 1);
    coord n = offset_normal (point, ibmf, xindex, yindex);
   
    double leftWeight = 0, rightWeight = 0, bottomWeight = 0, topWeight = 0;

    // sum of x contributions
    int leftIndex = xindex - 1, rightIndex = xindex + 1;
    if (ibm[leftIndex,yindex] > 0.5) {
        leftWeight = sq(n.x) * ibmf.x[xindex,yindex];
    }
    else if (ibm[rightIndex,yindex] > 0.5) { // should be else if? or separate if?
        rightWeight = sq(n.x) * ibmf.x[rightIndex,yindex];
    }

    // sum of y contributions
    int bottomIndex = yindex - 1, topIndex = yindex + 1;
    if (ibm[xindex,bottomIndex] > 0.5) {
        bottomWeight = sq(n.y) * ibmf.y[xindex,yindex];
    }
    else if (ibm[xindex,topIndex] > 0.5) { // should be else if? or separate if?
        topWeight = sq(n.y) * ibmf.y[xindex,topIndex];
    }

    // calculate flux of entire ghost cell
    coord nOutward, midPoint;
    double area = ibm_geometry (point, &midPoint, &nOutward);

    double veloFlux = -uf.x[xindex,yindex] * ibmf.x[xindex,yindex] +
                       uf.x[rightIndex,yindex] * ibmf.x[rightIndex,yindex] +
                      -uf.y[xindex,yindex] * ibmf.y[xindex,yindex] +
                       uf.y[xindex,topIndex] * ibmf.y[xindex,topIndex] -
                       uibm_x(midPoint.x,midPoint.y,midPoint.z) * nOutward.x * area - 
                       uibm_y(midPoint.x,midPoint.y,midPoint.z) * nOutward.y * area;
    // veloFlux /= Delta;

    double weightSum = leftWeight + rightWeight + bottomWeight + topWeight + SEPS;

    double xFlux = veloFlux * ((leftWeight + rightWeight) / weightSum);
    double yFlux = veloFlux * ((topWeight + bottomWeight) / weightSum);
    coord fluxes = {xFlux, yFlux};

    if (ibm[] <= 0.5) {
        foreach_dimension() {
            fluxes.x *= -1;
        }
    }
#if 0
    fprintf (stderr, "%g %g n.x=%g n.y=%g ws=%g f.x=%g f.y=%g\n",
                       x, y, n.x, n.y, weightSum, fluxes.x, fluxes.y);
#endif
    return fluxes;
}

foreach_dimension()
double virtual_merge_x (Point point, scalar ibm, face vector ibmf, face vector uf)
{
    if (ibm[] <= 0) {
        return 0;
    }
    int index = borders_ghost_x(point, ibm);

    // nothing special to be done for cells that don't border ghost cells
    if (!index && (ibm[] > 0.5 || ibm[] <= 0)) {
        return 0;
    }

    coord mergedFlux = ghost_fluxes (point, ibm, ibmf, uf);
    
    return mergedFlux.x;
}
#endif



/*
The next few functions are taken from embed to calculate interfacial force.
*/

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar ibm,
					   coord n, coord p, double bc,
					   double * coef)
{
  //foreach_dimension()
  //  n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !ibmf.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (ibmf.x[i + (i < 0),j] && ibmf.y[i,j] && ibmf.y[i,j+1] &&
	  ibm[i,j-1] && ibm[i,j] && ibm[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = ibmf.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!ibmf.y[i,j,k+m] || !ibmf.y[i,j+1,k+m] ||
	    !ibmf.z[i,j+m,k] || !ibmf.z[i,j+m,k+1] ||
	    !ibm[i,j+m,k-1] || !ibm[i,j+m,k] || !ibm[i,j+m,k+1])
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

  //fprintf(stderr, "(%g,%g) d0=%g d1=%g bc=%g v0=%g v1=%g ibm=%g\n", x, y, d[0], d[1], bc, v[0], v[1], ibm[]);

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

double dirichlet_gradient (Point point, scalar s, scalar ibm,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, ibm, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, ibm, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, ibm, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, ibm, n, p, bc, coef);
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
            dudn.x = dirichlet_gradient (point, u.x, ibm, n, p, vb, &val);
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


/*
ibm_force calculates the pressure force, Fp, and viscous force, Fmu, acting on the
immersed boundary.
*/

void ibm_force (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};
    foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        // if cell contains boundary intercept
        if (ibm[] > 0. && ibm[] < 1.) {
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
            double boundaryPressure = extrapolate_scalar (point, ibm, midPoint, n, p);
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


/*
This function fills cf with the skin friction coefficient for each interfacial cell.
// TODO: this assume no-slip condition; make it work for all B.C!
*/

double skin_friction (vector u, face vector mu, scalar cf)
{
    double cftotal = 0;
    foreach (reduction(+:cftotal)) {
        if (ibm[] > 0 && ibm[] < 1) {
            coord n, b;
            double area = ibm_geometry (point, &b, &n);
            area *= pow (Delta, dimension - 1);
            
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
            //coord tt = {-n.y, n.x}; // tangent vector
            //coord dudt = ibm_gradient (point, u, b, tt);
            foreach_dimension()
                tau.x -= mua* (dudn.x * (sq(n.x) + 1.) + 
                               dudn.y * -n.x * -n.y);
            cf[] = distance(tau.x, tau.y);
#else
            foreach_dimension()
                tau.x -= mua * (dudn.x * (sq(n.x) + 1.) + 
                                dudn.y * -n.x * -n.y +
                                dudn.z * -n.x * -n.z);
            cf[] = distance3D(tau.x, tau.y, tau.z);
#endif
            cftotal += cf[];
        }
        else
            cf[] = 0;
    }
    return cftotal;
}

/*
bilinear_ibm is another function taken from embed. If the parent cell is completely
inside the solib boundary, then simple injection is performed (i.e., the value of 
the child is set to that of the parent).

This is used in the multigrid solver, and is found to significantly improve 
convergence of the pressure solver.
*/

#if MULTIGRID
static inline double bilinear_ibm (Point point, scalar s)
{
    if (!coarse(ibm) || !coarse(ibm,child.x)) {
        return coarse(s);
    }
    #if dimension >= 2
    if (!coarse(ibm,0,child.y) || !coarse(ibm,child.x,child.y)) {
        return coarse(s);
    }
    #endif
    #if dimension >= 3
    if (!coarse(ibm,0,0,child.z) || !coarse(ibm,child.x,0,child.z) ||
        !coarse(ibm,0,child.y,child.z) ||
        !coarse(ibm,child.x,child.y,child.z)) {
        return coarse(s);  
    }
    #endif
 
    return bilinear (point, s);
}

#define bilinear(point, s) bilinear_ibm(point, s)
#endif // MULTIGRID



static void gradients_ibm (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
      if (s.gradient)
	foreach_dimension() {
#if IBM
      if (!ibmf.x[] || !ibmf.x[1])
        v.x[] = 0.;
      else
#endif
	    v.x[] = s.gradient (s[-1], s[], s[1])/Delta;
	}
      else // centered
	foreach_dimension() {
#if IBM
      if (!ibmf.x[] || !ibmf.x[1])
        v.x[] = 0.;
      else
#endif
	    v.x[] = (s[1] - s[-1])/(2.*Delta);
	}
    }
  }
}


static inline double vertex_average (Point point, scalar s)
{
#if dimension == 2
    return (4.*s[] + 
	        2.*(s[0,1] + s[0,-1] + s[1,0] + s[-1,0]) +
    	    s[-1,-1] + s[1,-1] + s[1,1] + s[-1,1])/16.;
#else
    return (8.*s[] +
	        4.*(s[-1] + s[1] + s[0,1] + s[0,-1] + s[0,0,1] + s[0,0,-1]) +
	        2.*(s[-1,1] + s[-1,0,1] + s[-1,0,-1] + s[-1,-1] + 
    		s[0,1,1] + s[0,1,-1] + s[0,-1,1] + s[0,-1,-1] +
	    	s[1,1] + s[1,0,1] + s[1,-1] + s[1,0,-1]) +
	        s[1,-1,1] + s[-1,1,1] + s[-1,1,-1] + s[1,1,1] +
    	    s[1,1,-1] + s[-1,-1,-1] + s[1,-1,-1] + s[-1,-1,1])/64.;
#endif
}


/*
local_to_global fills ax, ay, and az with the global coordinates of a point, p,
given in a local coordinate system.
*/

int local_to_global (Point point, coord p, double* ax, double* ay, double* az)
{
    *ax = x + p.x * Delta;
    *ay=  y + p.y * Delta;
    *az = z + p.z * Delta;

    return 0;
}


/*
copy_coord is used to fill three variables with the coresponding components of p.
This is used to avoid having the .x or _x indicies being automatically changed
within a foreach_dimension()
*/

int copy_coord (coord p, double* ax, double* ay, double* az)
{
    *ax = p.x;
    *ay = p.y;
    *az = p.z;

    return 1;
}


foreach_dimension()
double ibm_flux_x (Point point, scalar s, face vector mu, double * val)
{
    *val = 0.;
    if (ibm[] >= 1. || ibm[] <= 0.)
        return 0.;


    coord n = facet_normal (point, ibm, ibmf), mp;
    double alpha = plane_alpha (ibm[], n);
    double area = plane_area_center (n, alpha, &mp);
    normalize (&n);
    foreach_dimension()
        n.x *= -1;

    //ibm_geometry(point, &n, &mp); // why doesn't this work?

    double mpx, mpy, mpz;
    local_to_global(point, mp, &mpx, &mpy, &mpz);

    //double bc = uibm_x(mpx, mpy, mpz);
    bool dirichlet = false;
    double bc = s.boundary[immersed] (point, point, s, &dirichlet);
    
    double coef = 0.;
    //fprintf(stderr, "\n| ibm_flux (%g, %g) ibm=%g s=%g n.x=%g n.y=%g mp.x=%g mp.y=%g bc=%g\n",
    //                    x, y, ibm[], s[], n.x, n.y, mp.x, mp.y, bc);
    double grad = dirichlet_gradient (point, s, ibm, n, mp, bc, &coef);
    double mua = 0., fa = 0.;
    foreach_dimension() {
        mua += mu.x[] + mu.x[1];
        fa += ibmf.x[] + ibmf.x[1];
    }

    *val = - mua/(fa + SEPS)*grad*area/Delta;
    return - mua/(fa + SEPS)*coef*area/Delta;
}

#if 1 // testing some different interpolation functions like embed
#define ibm_avg(a,i,j,k)							\
  ((a[i,j,k]*(1.5 + ibm[i,j,k]) + a[i-1,j,k]*(1.5 + ibm[i-1,j,k]))/	\
   (ibm[i,j,k] + ibm[i-1,j,k] + 3.))

#if dimension == 2

#define face_condition(ibmf, ibm)						\
  (ibmf.x[i,j] > 0.5 && ibmf.y[i,j + (j < 0)] && ibmf.y[i-1,j + (j < 0)] &&	\
   ibm[i,j] && ibm[i-1,j])

foreach_dimension()
static inline double ibm_face_gradient_x (Point point, scalar a, int i)
{
  int j = sign(ibmf.x[i,1] - ibmf.x[i,-1]);
  assert (ibm[i] && ibm[i-1]);
  if (face_condition (ibmf, ibm))
    return ((1. + ibmf.x[i])*(a[i] - a[i-1]) +
	    (1. - ibmf.x[i])*(a[i,j] - a[i-1,j]))/(2.*Delta);
  return (a[i] - a[i-1])/Delta;
}

foreach_dimension()
static inline double ibm_face_value_x (Point point, scalar a, int i)
{
  int j = sign(ibmf.x[i,1] - ibmf.x[i,-1]);
  return face_condition (ibmf, ibm) ?
    ((1. + ibmf.x[i])*ibm_avg(a,i,0,0) + (1. - ibmf.x[i])*ibm_avg(a,i,j,0))/2. :
    ibm_avg(a,i,0,0);
}

#else // dimension == 3

foreach_dimension()
static inline coord embed_face_barycentre_z (Point point, int i)
{
  // Young's normal calculation
  coord n1 = {0};
  double nn = 0.;
  scalar f = ibmf.z;
  foreach_dimension(2) {
    n1.x = (f[-1,-1,i] + 2.*f[-1,0,i] + f[-1,1,i] -
	    f[+1,-1,i] - 2.*f[+1,0,i] - f[+1,1,i]);
    nn += fabs(n1.x);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  foreach_dimension(2)
    n1.x /= nn;
  // Position `p` of the face barycentre
  coord n, p1, p;
  ((double *)&n)[0] = n1.x, ((double *)&n)[1] = n1.y;
  double alpha = line_alpha (f[0,0,i], n);
  line_center (n, alpha, f[0,0,i], &p1);
  p.x = ((double *)&p1)[0], p.y = ((double *)&p1)[1], p.z = 0.;
  return p;
}

#define face_condition(ibmf, ibm)						\
  (ibmf.x[i,j,k] > 0.5 && (ibmf.x[i,j,0] > 0.5 || ibmf.x[i,0,k] > 0.5) &&	\
   ibmf.y[i,j + (j < 0),0] && ibmf.y[i-1,j + (j < 0),0] &&			\
   ibmf.y[i,j + (j < 0),k] && ibmf.y[i-1,j + (j < 0),k] &&			\
   ibmf.z[i,0,k + (k < 0)] && ibmf.z[i-1,0,k + (k < 0)] &&			\
   ibmf.z[i,j,k + (k < 0)] && ibmf.z[i-1,j,k + (k < 0)] &&			\
   ibm[i-1,j,0] && ibm[i-1,0,k] && ibm[i-1,j,k] &&				\
   ibm[i,j,0] && ibm[i,0,k] && ibm[i,j,k])

foreach_dimension()
static inline double ibm_face_gradient_x (Point point, scalar a, int i)
{
  assert (ibm[i] && ibm[i-1]);
  coord p = embed_face_barycentre_x (point, i);
  // Bilinear interpolation of the gradient (see Fig. 1 of Schwartz et al., 2006)
  int j = sign(p.y), k = sign(p.z);
  if (face_condition(ibmf, ibm)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return (((a[i,0,0] - a[i-1,0,0])*(1. - p.y) +
	     (a[i,j,0] - a[i-1,j,0])*p.y)*(1. - p.z) + 
	    ((a[i,0,k] - a[i-1,0,k])*(1. - p.y) +
	     (a[i,j,k] - a[i-1,j,k])*p.y)*p.z)/Delta;
  }
  return (a[i] - a[i-1])/Delta;
}

foreach_dimension()
static inline double ibm_face_value_x (Point point, scalar a, int i)
{
  coord p = embed_face_barycentre_x (point, i);
  // Bilinear interpolation
  int j = sign(p.y), k = sign(p.z);
  if (face_condition(ibmf, ibm)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return ((ibm_avg(a,i,0,0)*(1. - p.y) + ibm_avg(a,i,j,0)*p.y)*(1. - p.z) + 
	    (ibm_avg(a,i,0,k)*(1. - p.y) + ibm_avg(a,i,j,k)*p.y)*p.z);
  }
  return ibm_avg(a,i,0,0);
}
#endif // dimension == 3

#if 0 // TODO: try turning third on/off and include it in macros
attribute {
    bool third;
}
#endif

#if 0
#undef face_gradient_x
#define face_gradient_x(a,i)					\
  (ibmf.x[i] < 1. && ibmf.x[i] > 0. ?			\
   ibm_face_gradient_x (point, a, i) :			\
   (a[i] - a[i-1])/Delta)

#undef face_gradient_y
#define face_gradient_y(a,i)					\
  (ibmf.y[0,i] < 1. && ibmf.y[0,i] > 0. ?		\
   ibm_face_gradient_y (point, a, i) :			\
   (a[0,i] - a[0,i-1])/Delta)

#undef face_gradient_z
#define face_gradient_z(a,i)					\
  (ibmf.z[0,0,i] < 1. && ibmf.z[0,0,i] > 0. ?		\
   embed_face_gradient_z (point, a, i) :			\
   (a[0,0,i] - a[0,0,i-1])/Delta)
#endif
#undef face_value
#define face_value(a,i)							\
  (true && ibmf.x[i] < 1. && ibmf.x[i] > 0. ?				\
   ibm_face_value_x (point, a, i) :					\
   ibm_avg(a,i,0,0))

#undef center_gradient
#define center_gradient(a) (ibmf.x[] && ibmf.x[1] ? (a[1] - a[-1])/(2.*Delta) : \
			    ibmf.x[1] ? (a[1] - a[])/Delta :		    \
			    ibmf.x[]  ? (a[] - a[-1])/Delta : 0.)
#endif


/*
The metric event is used to set the metric fields, fm and cm, to the ibmFaces and
ibmCells field, respectively. It is also used to specifiy the prolongation and
refinement operations for each field (the functions of which are defined in ibm-tree.h.

The definition for ibmFaces and ibmCells are as follows:

ibmCells = 0 if ghost or solid cell and 1 if fluid cell.
ibmFaces = 0 if it borders two solid cells and 1 if it borders two fluid cells
           or 1 fluid cell and 1 solid/ghost cell.
*/
#if TREE
#include "ibm-tree.h"
#endif
#if 1
event metric (i = 0)
{
    if (is_constant (fm.x)) {
        foreach_dimension()
            assert (constant (fm.x) == 1.);
        fm = ibmFaces;
    }
    foreach_face() {
        ibmFaces.x[] = 1.;
        ibmf.x[] = 1;
    }
    if (is_constant (cm)) {
        assert (constant (cm) == 1.);
        cm = ibmCells;
    }
    foreach() {
        ibmCells[] = 1.;
        ibm[] = 1.;
        ibm0[] = 1.;
    }

#if TREE
    // set prolongation and refining functions
    //ibmCells.refine = ibm_fraction_refine;
    //ibmCells.prolongation = fraction_refine;

    ibmCells.refine = fraction_refine_metric;
    ibmCells.prolongation = fraction_refine_metric;

    // THIS DOSEN'T WORK WITH AMR, EVEN IN SERIAL
    //ibmCells.restriction = restriction_cell_metric;

    ibm0.refine = ibm.refine = ibm_fraction_refine;
    ibm0.prolongation = ibm.prolongation = fraction_refine;

    foreach_dimension() {
        ibmFaces.x.prolongation = refine_metric_injection_x;
        ibmf.x.prolongation = ibm_face_fraction_refine_x;
        ibmFaces.x.restriction = restriction_face_metric; // is this really necessary?
    }
#endif
    restriction ({ibm, ibmf, ibmFaces, ibmCells});

    boundary(all);
}
#endif
