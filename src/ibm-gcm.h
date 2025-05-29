#undef dv
#define dv() (pow(Delta,dimension)*ibm[])

#include "fractions.h"

#define BGHOSTS 2
#define IBM 1
#define LIMIT 1.e10
#define INT_TOL 1e-6    // tolerance used for volume fraction fields (interface tolerance)

#undef SEPS
#define SEPS 1e-30

scalar ibm[];
scalar ibm0[];          // solid volume fraction field of previous timestep
face vector ibmf[];
face vector ibmf0[];

// metric fields
scalar ibmCells[];
face vector ibmFaces[];

// coord vc = {0,0,0};     // imposed velocity boundary condition (depreciated)


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


#define distance(a,b) sqrt(sq(a) + sq(b))
#define distance3D(a,b,c) sqrt(sq(a) + sq(b) + sq(c))
#define on_interface(a) (a[] > 0+INT_TOL && a[] < 1-INT_TOL)
#define is_mostly_solid(a, i) (a[i] > 0+INT_TOL && a[i] <= 0.5)
#define is_fresh_cell(a0, a) (a0[] <= 0.5 && a[] > 0.5)


/*
These functions and macros below are used to mimic Basilisk's way of specifying
boundary condtions, e.g., u.t[embed] = dirichlet(0). For now, we only allow a
dirichlet b.c for velocity.

TODO: allow neumann b.c (navier-slip b.c?)
TODO: make other macros for p, pf, f, etc.
TODO: use new macro feature
*/

static inline double uibm_x (double x, double y, double z);
static inline double uibm_y (double x, double y, double z);
#if dimension == 3
static inline double uibm_z (double x, double y, double z);
#endif

#define u_x_ibm_dirichlet(expr) \
    static inline double uibm_x (double x, double y, double z) {return expr;} \

#define u_y_ibm_dirichlet(expr) \
    static inline double uibm_y (double x, double y, double z) {return expr;} \

#define u_z_ibm_dirichlet(expr) \
    static inline double uibm_z (double x, double y, double z) {return expr;} \


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
   return fluid_neighbor(point, ibm) && ibm[] <= 0.5 && match_level(point, ibm);
}


/*
centroid_point returns the area of the interfrace fragment in a give cell. It
takes in the volume fraction field ibm and fills midPoint with the interfacial 
centroid in the GLOBAL coordinate system.

Note here n is the inward facing normal normalized so |n.x| + |n.y| + |n.z| = 1
*/

double centroid_point (Point point, scalar ibm, coord * midPoint, coord * n)
{
    coord cellCenter = {x, y, z};
    *n = facet_normal (point, ibm, ibmf);
    double alpha = plane_alpha (ibm[], *n);
    double area = plane_area_center (*n, alpha, midPoint);

    foreach_dimension()
        midPoint->x = cellCenter.x + midPoint->x*Delta;
    return area;
}


/*
The function below fills frag with the normal vector n, alpha, and the volume fraction of the
cell that is closest to the ghost cell. It also returns the coordinates of the fragment's midpoint 
and fills fluidCell with the cell center coordinates of the closest fluid cell.

Note: we make a crude approximation that the closest interfacial point from the surrounding
cells to the ghost cell is which ever interfacial mid point/centroid is closest. In practice, it has worked
adequately, but this can be improved.

TODO: Clean up and streamline function.
*/

coord closest_interface (Point point, vector midPoints, scalar ibm, 
                         vector normals, fragment * frag, coord * fluidCell)
{
    fragment temp_frag;
    coord temp_midPoint, temp_fluidCell = {0,0};
    coord n;
    double min_distance = 1e6;

     for(int i = -1; i <= 1; i++) {
        double dx = midPoints.x[i] - x;
        double dy = midPoints.y[i] - y;
        double dz = midPoints.z[i] - z;
        if (midPoints.x[i] && distance3D(dx, dy, dz) < min_distance) {
            temp_midPoint.x = midPoints.x[i];
            temp_midPoint.y = midPoints.y[i];
            temp_midPoint.z = midPoints.z[i];

            n.x = normals.x[i]; n.y = normals.y[i]; n.z = normals.z[i];

            fill_fragment (ibm[i], n, &temp_frag);
            temp_fluidCell.x = i*Delta + x;
            temp_fluidCell.y = y;
            temp_fluidCell.z = z;
            min_distance = distance3D(dx, dy, dz);
        }
     }

     for(int j = -1; j <= 1; j++) {
        double dx = midPoints.x[0,j] - x;
        double dy = midPoints.y[0,j] - y;
        double dz = midPoints.z[0,j] - z;
        if (midPoints.x[0,j] && distance3D(dx, dy, dz) < min_distance) {
            temp_midPoint.x = midPoints.x[0,j];
            temp_midPoint.y = midPoints.y[0,j];
            temp_midPoint.z = midPoints.z[0,j];

            n.x = normals.x[0,j]; n.y = normals.y[0,j]; n.z = normals.z[0,j];

            fill_fragment (ibm[0,j], n, &temp_frag);
            temp_fluidCell.x = x;
            temp_fluidCell.y = j*Delta + y;
            temp_fluidCell.z = z;
            min_distance = distance3D(dx, dy, dz);
        }
     }
#if dimension == 3
     for(int k = -1; k <= 1; k++) {
        double dx = midPoints.x[0,0,k] - x;
        double dy = midPoints.y[0,0,k] - y;
        double dz = midPoints.z[0,0,k] - z;
        if (midPoints.x[0,0,k] && distance3D(dx, dy, dz) < min_distance) {
            temp_midPoint.x = midPoints.x[0,0,k];
            temp_midPoint.y = midPoints.y[0,0,k];
            temp_midPoint.z = midPoints.z[0,0,k];

            n.x = normals.x[0,0,k]; n.y = normals.y[0,0,k]; n.z = normals.z[0,0,k];

            fill_fragment (ibm[0,0,k], n, &temp_frag);
            temp_fluidCell.x = x;
            temp_fluidCell.y = y;
            temp_fluidCell.z = k*Delta + z;
            min_distance = distance3D(dx, dy, dz);
        }
     }
#endif
    *fluidCell = temp_fluidCell;
    *frag = temp_frag;

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
    normalize(&n);
    // double mag = fabs(n.x) + fabs(n.y) + fabs(n.z);

    double offset = 0;
    offset += n.x * -sign2(fluidCell.x - x);
    offset += n.y * -sign2(fluidCell.y - y);
    offset += n.z * -sign2(fluidCell.z - z);
    coord boundaryInt = {(-frag.alpha / mag - offset) * n.x,
                         (-frag.alpha / mag - offset) * n.y,
                         (-frag.alpha / mag - offset) * n.z};


#if 0
    if (is_ghost_cell(point, ibm))
        fprintf (stderr, "\n || bi.x=%g bi.y=%g sum=%g\n",
                                boundaryInt.x, boundaryInt.y, offset);
#endif

    foreach_dimension()
        boundaryInt.x = ghostCell.x + boundaryInt.x*Delta;

#if 0
    if (is_ghost_cell(point, ibm)) {
        fprintf (stderr, "|| %g %g bi.x=%g bi.y=%g fluid.x=%g fluid.y=%g\n", 
                              x, y, boundaryInt.x, boundaryInt.y, fluidCell.x, fluidCell.y);
        fprintf (stderr, "|| n.x=%g n.y=%g alpha=%g mag=%g\n",
                             frag.n.x, frag.n.y, frag.alpha, mag);
    }
#endif
    return boundaryInt;
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

TODO: Signals when cell is on MPI and/or openMP boundary too?
TODO: 3D implementation
*/

bool borders_boundary (Point point, int * useri = NULL, int * userj = NULL)
{
#if TREE
    // Look at directly adjacent neighbors (4 in 2D)
    for (int d = 0; d < dimension; d++) {
	    for (int k = -1; k <= 1; k += 2) {
            int i = 0, j = 0;
	        if (d == 0) {
                i = k; 
            }
            else if (d == 1) {
                j = k; 
            }
            // check to see if neighboring cell is inside boundary
            if (neighbor(-i,-j,0).pid < 0) {

                if (useri) *useri = -i;
                if (userj) *userj = -j;

                return true;
            }
        }
    }
#endif
    return false;
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
*/

void fluid_only (Point point, int xx, int yy, int zz, int i, int j, int k, 
                 coord * pTemp, coord * velocity, vector midPoints,
                 int bOffset_X = 0, int bOffset_Y = 0, int bOffset_Z = 0)
{
    int off_x = xx + i, off_y = yy + j, off_z = zz + k;
    if (ibm[off_x,off_y,off_z] <= 0.5 && ibm[off_x,off_y,off_z] > 0.) {
        pTemp->x = midPoints.x[off_x,off_y,off_z];
        pTemp->y = midPoints.y[off_x,off_y,off_z];
        pTemp->z = midPoints.z[off_x,off_y,off_z];

        double mpx = pTemp->x, mpy = pTemp->y, mpz = pTemp->z;
        foreach_dimension() {
            if (bOffset_X == off_x) {
                pTemp->x += bOffset_X * Delta;
            }
            velocity->x = uibm_x(mpx, mpy, mpz);
        }
    }
}


/*
The function below uses interpolation to find the velocity at the image point and
returns it given a vector field, u, the coordinates of the image point, and a field
containing all interfacial midpoints.

TODO: Streamline and clean-up code?
TODO: Extend to handle higher-order interpolation schemes, i.e. larger matrices.
*/

coord image_velocity (Point point, vector u, coord imagePoint, vector midPoints)
{
    
    int boundaryOffsetX = 0, boundaryOffsetY = 0;
    bool border = borders_boundary (point, &boundaryOffsetX, &boundaryOffsetY);
    
    int xOffset = 0, yOffset = 0, zOffset = 0;
    image_offsets (point, imagePoint, &xOffset, &yOffset, &zOffset);
    
    assert (abs(xOffset) <= 2 && abs(yOffset) <= 2 && abs(zOffset) <= 2);

    coord imageCell = {x + Delta * xOffset, y + Delta * yOffset, z + Delta * zOffset};
    
    int i = sign(imagePoint.x - imageCell.x);
    int j = sign(imagePoint.y - imageCell.y);
    int k = sign(imagePoint.z - imageCell.z);
   
    #if 0
    if (border) {
        fprintf (stderr, "WARNING: cell x=%g y=%g borders a boundary: %d %d %d %d\n", 
                x, y, point.i, point.j, boundaryOffsetX, boundaryOffsetY);
    }
    #endif

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
                boundaryOffsetX, boundaryOffsetY);

    fluid_only (point, xx, yy, zz, i, 0, 0, &p1, &velocity[1], midPoints, 
                boundaryOffsetX, boundaryOffsetY);

    fluid_only (point, xx, yy, zz, i, j, 0, &p2, &velocity[2], midPoints, 
                boundaryOffsetX, boundaryOffsetY);

    fluid_only (point, xx, yy, zz, 0, j, 0, &p3, &velocity[3], midPoints, 
                boundaryOffsetX, boundaryOffsetY);
#if dimension == 3
    fluid_only (point, xx, yy, zz, 0, 0, k, &p4, &velocity[4], midPoints, 
                boundaryOffsetX, boundaryOffsetY);

    fluid_only (point, xx, yy, zz, i, 0, k, &p5, &velocity[5], midPoints, 
                boundaryOffsetX, boundaryOffsetY);

    fluid_only (point, xx, yy, zz, i, j, k, &p6, &velocity[6], midPoints, 
                boundaryOffsetX, boundaryOffsetY);

    fluid_only (point, xx, yy, zz, 0, j, k, &p7, &velocity[7], midPoints, 
                boundaryOffsetX, boundaryOffsetY);
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
    return temp_velo;

}


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
double ibm_geometry (Point point, coord * p, coord * n, double * alphau = NULL)
{
    *n = facet_normal (point, ibm, ibmf);
    double alpha = plane_alpha (ibm[], *n);
    if (alphau != NULL)
        *alphau = alpha;
    double area = plane_area_center (*n, alpha, p);
    foreach_dimension()
        n->x *= -1;
    normalize (n);

    return area;
}


static inline
double ibm0_geometry (Point point, coord * p, coord * n)
{
    *n = facet_normal (point, ibm0, ibmf0);
    double alpha = plane_alpha (ibm0[], *n);
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

TODO: clean up and streamline code
TODO: 3D implementation
TODO: what if ghost cell has two fluid cell neighbors in one direction
      (left and right or top and bottom)?
*/

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
    coord cellCenter = {x,y,z}, midPoint, dudn;
    foreach_dimension() {
        midPoint.x = cellCenter.x + p.x*Delta;
     }
    double px = midPoint.x, py = midPoint.y, pz = midPoint.z;
    foreach_dimension() {
        double vb = uibm_x(px,py,pz);
        double val;
        dudn.x = dirichlet_gradient (point, u.x, ibm, n, p, vb, &val);
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
            area *= pow (Delta, dimension - 1);

            
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

    double bc = uibm_x(mpx, mpy, mpz);
    //double bc = 0;
    
    double coef = 0.;
    //fprintf(stderr, "\n| ibm_flux (%g, %g) ibm=%g s=%g n.x=%g n.y=%g mp.x=%g mp.y=%g bc=%g\n",
    //                    x, y, ibm[], s[], n.x, n.y, mp.x, mp.y, bc);
    double grad = dirichlet_gradient (point, s, ibm, n, mp, bc, &coef);
    double mua = 0., fa = 0.;
    foreach_dimension() {
        mua += mu.x[] + mu.x[1];
        fa += fm.x[] + fm.x[1];
        //fa += ibmf.x[] + ibmf.x[1];
    }

    *val = - mua/(fa + SEPS)*grad*area/Delta;
    return - mua/(fa + SEPS)*coef*area/Delta;
}




#if 0 // this seems to not have any major effect, despite its use in embed
#define face_condition(ibmf, ibm)						\
  (ibmf.x[i,j] > 0.5 && ibmf.y[i,j + (j < 0)] && ibmf.y[i-1,j + (j < 0)] &&	\
   ibm[i,j] && ibm[i-1,j])

foreach_dimension()
static inline double ibm_face_gradient_x (Point point, scalar a, int i)
{
  if (ibmf.x[i] < 1. && ibmf.x[i] > 0.) {
    int j = sign(ibmf.x[i,1] - ibmf.x[i,-1]);
    assert (ibm[i] && ibm[i-1]);
    if (face_condition (ibmf, ibm))
      return ((1. + ibmf.x[i])*(a[i] - a[i-1]) +
	          (1. - ibmf.x[i])*(a[i,j] - a[i-1,j]))/(2.*Delta);
  }
  else
    return (a[i] - a[i-1])/Delta;
}
#endif


/*
###### TWO PHASE FUNCTIONS ######

TODO: move these functions to a separate header file
*/

(const) scalar contact_angle;

/*
boundary_points is used to find the intersecting points of the fluid interface
within the area being advected.
*/

int boundary_points (coord nf, double alphaf, coord lhs, coord rhs, coord bp[2])
{
    int i = 0;
    // check on x faces first
    double dx = rhs.x - lhs.x;
    if (fabs(dx) < 1e-15)
        return 0;

    if (fabs(nf.x) < 1e-15 && fabs(nf.y) < 1e-15)
        return 0;

    for (double xint = rhs.x; xint >= lhs.x; xint -= dx) {
        double yint = (alphaf - nf.x*xint)/(nf.y+SEPS);
        if (fabs(yint) <= 0.5) {
            bp[i].x = xint;
            bp[i].y = yint;
            ++i;
        }
    }

    // then check y faces
    for (double yint = -0.5; yint <= 0.5; yint += 1.) {
        double xint = (alphaf - nf.y*yint)/(nf.x+SEPS);
        if (xint <= rhs.x && xint >= lhs.x) {
            bp[i].x = xint;
            bp[i].y = yint;
            ++i;
        }
    }

    return i;
}


/*
this function checks to see if the two interfaces intersect each other. If they
do, and it is within the bounds of the region, the coordinates is stored in pi.
It also returns 1 or 0 based on if the lines intersect each other or not.
*/

int interface_intersect (coord nf, double alphaf, coord ns, double alphas,
                         coord lhs, coord rhs, coord * pint)
{
    coord pt;
    pt.x = (alphas/(ns.y + SEPS) - alphaf/(nf.y + SEPS)) /
                  ((ns.x/(ns.y + SEPS)) - (nf.x/(nf.y + SEPS)) + SEPS);

    pt.y = (alphaf/(nf.y + SEPS)) - (nf.x*pt.x)/(nf.y + SEPS);

    foreach_dimension() {
        if (pt.x < lhs.x || pt.x > rhs.x) {
            return 0;
        }
    }

    *pint = pt;
    return 1;
}


/*
this function checks to see if a given point, pc, is inside the region
containing only real fluid and not inside the immersed boundary.

Note: ns is the inward pointing normal for the solid boundary while
      nf is the outward pointing normal for the fluid boundary.

returns
    + if inside, - if outisde, and 0 if the points is on the interface
*/

double region_check (double vol, coord pc, coord nf, double alphaf, coord ns, double alphas)
{
    double fluid = alphaf - nf.x*pc.x - nf.y*pc.y;
    double solid = alphas - ns.x*pc.x - ns.y*pc.y;

    return difference(fluid,-solid); //- because solid has inward pointing normal
    //return intersection(fluid,solid); // these options are equivalent
    //return min(fluid,solid);
}


/*
This function uses cross products to determine the orientation of two points w.r.t
a given center coordinate.

is_begin will return
    true is a behind b in a clockwise order,
    false if a is in front of b, i.e. a is already in clockwise order with b.

in degenerate cases where a and b lay along the same line (a x b = 0) we return
true if a is further from the center than b.
*/

bool is_behind (coord a, coord b, coord center)
{
    // calculate cross product of a and b w.r.t the center
    double det = (a.x - center.x)*(b.y - center.y) - (a.y - center.y)*(b.x - center.x);

    if (det > 0)
        return false;
    if (det < 0)
        return true;

    // a and b lay along the same line, so det = 0. Instead, we check to see
    // if a is further from the center than b.
    double dista = sq(a.x - center.x) + sq(a.y - center.y);
    double distb = sq(b.x - center.x) + sq(b.y - center.y);

    return dista > distb;
}


/*
sort_clockwise sorts a list of coordinates, provided in cf w/nump points, in
clockwise order (or counter-clockwise if y-advection).
*/

void sort_clockwise (int nump, coord cf[nump])
{
    double xsum = 0, ysum = 0;
    for (int i = 0; i < nump; ++i) {
        xsum += cf[i].x;
        ysum += cf[i].y;
    }
    coord pc = {xsum/nump, ysum/nump}; // center coordinate is average of all points
    
    bool sorted = false;
    int nitrmax = 20, nitr = 0;
    while (!sorted && nitr < nitrmax) {
        sorted = true;
        for (int i = 0; i < nump; ++i) {
            int next = i + 1 < nump? i + 1: 0; // we loop around the list of coords
            if (is_behind(cf[i],cf[next],pc)) {
                coord tempcf = cf[i];
                cf[i] = cf[next];
                cf[next] = tempcf;
                sorted = false;
            }
        }
        nitr++;
    }
}


/*
polygon_area calculates the area enclosed by a list of points (in cw or ccw order)
using the shoelace formula.
*/

double polygon_area (int nump, coord cf[nump])
{
    double area = 0;
    bool finished = false;
    for (int i = 0; i < nump; ++i) {
        int next = i + 1 < nump? i + 1: 0; // to close the shape
        area += cf[i].x*cf[next].y - cf[next].x*cf[i].y;
    }

    return fabs(area)/2.;
}


/*
line_intersect returns the x (resp. y) coordinate along a line given y (resp. x).
The return value is in the cell's local coordinate system, i.e. [-0.5,-0.5]x[0.5,0.5].
*/

double line_intersect (double alpha, coord n, double x = HUGE, double y = HUGE)
{
    double inter = 0;
    if (x == HUGE && y != HUGE)
        inter = alpha/n.x - (n.y/(n.x+SEPS))*y;
    else
        inter = alpha/n.y - (n.x/(n.y+SEPS))*x;
    return inter;
}

// adjust x coordinate (rhsx) to conserve area
// lhsb = left-hand-side bottom, rhst = right-hand-side top
double fit_area (coord ns, double alphas, coord lhsb, coord rhst)
{
    double oldx = rhst.x;
    // assuming liquid is on top of solid interface
    coord lhst = {lhsb.x, rhst.y}, rhsb = {rhst.x, lhsb.y};

    // change the verticies if the solid interface is on top of the liquid interface
    if (ns.y > 0) {
        rhst.y = lhsb.y;
        lhsb.y = -0.5;
        lhst.y = rhst.y;
    }

    double error = HUGE, tolerance = 1e-9;
    int itr = 0, maxitr = 40;

    coord rect0[4] = {lhsb,rhsb,rhst,lhst};
    double area0 = polygon_area (4, rect0);
    double areaReal = area0 * rectangle_fraction (ns, alphas, lhsb, rhst);
    double error0 = area0 - areaReal;

    double minx = -0.5, maxx = 0, newx; // maxx could probably be 0 instead due to CFL=0.5 (was 0.5)
    while (fabs(error) > tolerance && itr < maxitr) {

        newx = (maxx + minx)/2.; // find mid section

        coord lhsb2 = {lhsb.x,-0.5}, rhsb2 = {newx, -0.5}, rhst2 = {newx, 0.5}, lhst2 = {lhst.x,0.5};
        coord rect[4] = {lhsb2,rhsb2,rhst2,lhst2};

        double areaTotal = polygon_area (4, rect);
        areaReal = areaTotal * rectangle_fraction (ns, alphas, lhsb2, rhst2);

        error = area0 - areaReal;
        if (error > 0)
            minx = newx;
        else
            maxx = newx;
        //fprintf(stderr, "AREA SOLVER: %d error=%g area0=%g arear=%g minx=%g maxx=%g\n",
        //                 itr, error, area0, areaReal, minx, maxx);
        itr++;
    }
    fprintf(stderr, "AREA SOLVER done: oldx=%g newx=%g\n", oldx, newx);
    if (itr == maxitr)
        fprintf(stderr, "WARNING: area root solver didn't converge after %d iteration wtih an error of %g\n",
                        itr, error);
    return newx;
}


/*
immersed_area calculates the area of the real fluid part of a cell.
*/

double immersed_area (double c, coord nf, double alphaf, coord ns, double alphas, 
                      coord lhs, coord rhs, int print = 0)
{
    if (lhs.x == rhs.x)
        return 0;

    // 1. calculate the area of the real region within the advected region
    coord lhst = {lhs.x,0.5}, rhsb = {rhs.x,-0.5}; 
    coord rect[4] = {lhs,rhsb,rhs,lhst};
    double areaTotal = polygon_area (4, rect); // total area being considered for advection (uf*dt*h)
    double areaLiquid = rectangle_fraction (ns, alphas, lhs, rhs); // volume fraction that isn't solid

    double cvy = line_intersect (alphas, ns, x = -0.5), rhsx = 0; // intercept of solid interface on left face
    #if 1
    if (rhs.x < 0.5 && areaLiquid < 1) {
        cvy = clamp(cvy, -0.5, 0.5);
        rhsx = fit_area (ns, alphas, (coord){lhs.x, cvy}, (coord){rhs.x, rhs.y});
        rhs.x = rhsx;
        rhsb.x = rhsx;
        coord rect0[4] = {lhs,rhsb,rhs,lhst};
        areaTotal = polygon_area (4, rect0);
        areaLiquid = rectangle_fraction (ns, alphas, lhs, rhs);
    }
    #endif
    double areaAdv = areaTotal*areaLiquid;

    // 2. find the intersection points, pf & ps, of the fluid and solid interface
    //    with the enclosed region
    coord pf[2], ps[2];
    for (int i = 0; i < 2; ++i)
        foreach_dimension() {
            pf[i].x = nodata;
            ps[i].x = nodata;
        }

    int numpf = boundary_points(nf, alphaf, lhs, rhs, pf);
    int numps = boundary_points(ns, alphas, lhs, rhs, ps);

    // 3. find the intersecting point of the two interfaces (if there is one)
    coord pint = {nodata,nodata,nodata};
    int numpi = interface_intersect (nf, alphaf, ns, alphas, lhs, rhs, &pint);


    // 4. find which points create the polygon defining the real fluid region
    //    9 possible points
    coord poly[9];
    poly[0] = pf[0], poly[1] = pf[1], poly[2] = ps[0], poly[3] = ps[1], poly[4] = pint;
    poly[5] = lhs, poly[6] = rhs, poly[7] = lhst, poly[8]= rhsb;

    int nump = 0; // # of real points
    for (int i = 0; i < 9; ++i) {
        double placement = region_check(c, poly[i], nf, alphaf, ns, alphas);
        if ((placement >= 0 || fabs(placement) < 1e-6) && poly[i].x != nodata) // should be < nodata?
            nump++;
        else 
            poly[i].x = nodata, poly[i].y = nodata;
    }

    if (print == 1) {
        fprintf(stderr, "||  pf1=(%g,%g) pf2(%g,%g) ps1=(%g,%g) ps2=(%g,%g) pint=(%g,%g)\n",
                        pf[0].x, pf[0].y, pf[1].x, pf[1].y, ps[0].x, ps[0].y, ps[1].x, 
                        ps[1].y, pint.x, pint.y);
   
        fprintf(stderr, "|| p0=(%g,%g) p1=(%g,%g) p2=(%g,%g) p3=(%g,%g) p4=(%g,%g)"
                        " p5=(%g,%g) p6=(%g,%g) p7=(%g,%g) p8=(%g,%g)\n", 
                        poly[0].x, poly[0].y, poly[1].x, poly[1].y,poly[2].x, poly[2].y, 
                        poly[3].x, poly[3].y, poly[4].x, poly[4].y, poly[5].x, poly[5].y, 
                        poly[6].x, poly[6].y, poly[7].x, poly[7].y, poly[8].x, poly[8].y);

        fprintf(stderr, "||  lhs=(%g,%g) rhs=(%g,%g) lhst=(%g,%g) rhsb=(%g,%g)\n",
                        lhs.x, lhs.y, rhs.x, rhs.y, lhst.x, lhst.y, rhsb.x, rhsb.y);
                        
        fprintf(stderr, "||  nf=(%g,%g) alphaf=%g ns=(%g,%g) alphas=%g c=%g %d %d %d %d\n",
                        nf.x, nf.y, alphaf, ns.x, ns.y, alphas, c, numpf, numps, numpi, nump);
    }

    if (nump == 0) // we don't do anything if the region has no interface fragments
        return 0;

    coord cf[nump]; // holds real points
    int count = 0;
    for (int i = 0; i < 9; ++i)
        if (poly[i].x != nodata && count < nump) {
            cf[count] = poly[i];
            count++;
        }

    if (print == 1) 
        for (int i = 0; i < nump; ++i) 
            fprintf(stderr, "cf[%d] = (%g,%g)\n", i, cf[i].x, cf[i].y);

    // 5. sort the real points in clockwise order
    sort_clockwise (nump, cf);

    // 6. use the shoelace formula to find the area
    double area = polygon_area (nump, cf);
    // coord rect[4] = {lhs,rhsb,rhs,lhst};
    // double areaTotal = polygon_area (4, rect);
    // double areaLiquid = rectangle_fraction (ns, alphas, lhs, rhs);

    #if 0
    if (fabs(cvy) < 0.5) {
        coord lhs2 = {lhs.x, cvy}, rhs2 = {rhsb.x, cvy};
        if (ns.y > 0) {
            coord rect[4] = {lhs, rhsb, rhs2, lhs2};
            areaAdv = polygon_area (4, rect);
            if (print == 1)
                for (int i = 0; i < 4; ++i)
                    fprintf(stderr, "(a) rect[%d] = (%g,%g)\n",
                             i, rect[i].x, rect[i].y);
        }
        else {
            coord rect[4] = {lhs2, rhs2, rhs, lhst};
            areaAdv = polygon_area (4, rect);
            if (print == 1)
                for (int i = 0; i < 4; ++i)
                    fprintf(stderr, "(b) rect[%d] = (%g,%g)\n",
                             i, rect[i].x, rect[i].y);
        }
    }
    else
        areaAdv = areaTotal*areaLiquid;
    #endif
    if (print == 1) {
        fprintf (stderr, "AFTER SORTING\n");
        for (int i = 0; i < nump; ++i) {
            fprintf(stderr, "cf[%d] = (%g,%g)\n",
                         i, cf[i].x, cf[i].y);
        }
        // get f[] w/o considering immersed boundary
        double f0 = rectangle_fraction(nf, alphaf, lhs, rhs);
        fprintf (stderr, "area=%0.15g  areaTotal=%0.15g areaLiquid=%0.15g "
                         "areaf=%0.15g f0=%0.15g areaAdv=%g cvy=%g vf=%g\n", 
                          area, areaTotal, areaLiquid, area/(areaTotal*areaLiquid), 
                          f0, areaAdv, cvy, area/areaAdv);
    }

    double vf = clamp(area/areaAdv, 0., 1.);
    //return vf >= 1 - INT_TOL? 1.: vf; // necessary for incline?
    return vf;
}

/*
This function calculates the fraction of a rectangle (defined by lhs and rhs)
which lies inside the liquid interface neglecting the portion inside of the 
immersed boundary.

lhs and rhs are the bottom left and top right (resp.) coordinates defining the
region being advected by the split VOF advection scheme (see sweep_x in vof.h)
*/

double immersed_fraction (double c, coord nf, double alphaf, coord ns, double alphas, 
                          coord lhs, coord rhs, int print = 0)
{
    return immersed_area(c, nf, alphaf, ns, alphas, lhs, rhs, print);
}



/*
immersed_line_alpha calculates the alpha value that conserves the volume of
real fluid, given in freal. We find the root of a function using the iterative
bisection method to obtain alpha.

The solver converges after about 20 iterations if tolerance = 1e-7

TODO: implement proper min and max solver for alpha
*/

#define PRINTA 1


typedef struct tripoint
{
    coord nf, ns;
    double alphaf, alphas;
    double f, fr, s;

} tripoint;


tripoint fill_tripoint (double fr, coord nf, double alphaf, coord ns, double alphas,
                        double f = 0, double s = 0)
{
    tripoint tcell = {nf, ns, alphaf, alphas, f, fr, s};
    return tcell;
}


int bisection_solver (double* a, double amin, double amax, tripoint tcell, 
                      double tolerance = 1e-9, int maxitr = 40)
{
    coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}; // bottom-left & top-right points, resp.
   
    double ibm0 = tcell.s == 0? rectangle_fraction (tcell.ns, tcell.alphas, lhs, rhs): tcell.s;

    double cmax = rectangle_fraction (tcell.nf, amax, lhs, rhs);
    double crmax = immersed_fraction (cmax, tcell.nf, amax, tcell.ns, tcell.alphas, lhs, rhs);
    double errmax = crmax - tcell.fr;
    //if (fabs(errmax) <= tolerance) {
    //    *a = tcell.alphaf;
    //    return 0;
    //}

    double alpha = 0, error = HUGE;
    int itr = 0;
    while (fabs(error) > tolerance && itr < maxitr) {
        alpha = (amin + amax)/2.;   // bisection
        double fa = rectangle_fraction (tcell.nf, alpha, lhs, rhs);
        double fcalc = immersed_fraction (fa, tcell.nf, alpha, tcell.ns, tcell.alphas, lhs, rhs);
        error = fcalc - tcell.fr;
        if (sign2(error) == sign2(errmax)) {
            amax = alpha;
            errmax = error;
          }
        else
            amin = alpha;
        itr++;

          #if PRINTA
           fprintf(stderr, "|| A.S: %d amin=%0.15g amax=%0.15g alpha=%0.15g" 
                           " fa=%0.15g fcalc=%0.15g freal=%g error=%g maxerror=%g\n",
                           itr, amin, amax, alpha, fa, fcalc, tcell.fr, error, errmax);
          #endif
    }
        
    if (itr == maxitr) 
        fprintf(stderr, "WARNING: alpha  solver does not converge after"
                        " maximum iteration (%d), error = %g\n", maxitr, error);

    *a = alpha;
    return itr;
}


int golden_search (double* a, double amin, double amid, double amax, tripoint tcell,
                   double tolerance = 1e-9, int maxitr = 40)
{
    coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}; // bottom-left & top-right points, resp.

    double ibm0 = tcell.s == 0? rectangle_fraction (tcell.ns, tcell.alphas, lhs, rhs): tcell.s;

    double alpha = 0, error = HUGE;
    int itr = 0;
    while (fabs(error) > tolerance && itr < maxitr) {
        double fmax = rectangle_fraction (tcell.nf, amax, lhs, rhs);
        double fmid = rectangle_fraction (tcell.nf, amid, lhs, rhs);
        double fmin = rectangle_fraction (tcell.nf, amin, lhs, rhs);

#if 0
        // if the cell has real fluid, we want to minimize f (think CA > 90)
        if (tcell.fr > 0) {
            if (fmax <= fmid && fmax <= fmin) alpha = fmax;
        }
        else { // the cell doesn't have real fluid, so we want to maximize f (think CA < 90)
        
        }
#endif
        if (fmax >= fmid && fmax >= fmin)
            alpha = amax;
        else if (fmid >= fmax && fmid >= fmin)
            alpha = amid;
        else if (fmin >= fmax && fmin >= fmid)
            alpha = amin;

        #if 0
        if (fmax >= 1-1e-7 && fmid < 1-1e-7)
            alpha = amid;
        else if (fmax >= 1-1e-6 && fmid >= 1-1e-6)
            alpha = amin;
        #else
         if (fmax >= 1 && fmid < 1)
            alpha = amid;
        else if (fmax >= 1 && fmid >= 1)
            alpha = amin;

        fprintf (stderr, "fmin=%g fmid=%g fmax=%g alpha=%g\n",
                          fmin, fmid, fmax, alpha);
        
        //if (rectangle_fraction (tcell.nf, alpha, lhs, rhs) < 1e-6)
        //    alpha = amin; // temp fix

        break;
    }
    if (itr == maxitr) 
    fprintf(stderr, "WARNING: alphaf (d) solver does not converge after"
                    " maximum iteration (%d), error = %g\n", maxitr, error);

    *a = alpha;
    return itr;
}


double immersed_line_alpha (Point point, coord nf, double alphaf, coord ns, double alphas,
                            double freal, double tolerance = 1e-9)
{
    if ((freal <= 1e-6 || freal >= 1-1e-6) && !nf.x && !nf.y)
        return alphaf;

    coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}; // bottom-left & top-right points, resp.

    double alphaMin = -1, alphaMax = 1; // are there better values?
 
    double error = HUGE, alpha = 0;

    double f0 = clamp(rectangle_fraction (nf, alphaf, lhs, rhs), 0., 1.);
    //if (f0 == 0. && alphaf == 0.)
    //    f0 = 1.;
    double ibm0 = rectangle_fraction (ns, alphas, lhs, rhs);
    freal = clamp(freal, 0., 1.);
    tripoint tcell = fill_tripoint (freal, nf, alphaf, ns, alphas, f0, ibm0);

    int maxitr = 40;

    // the cell is full or empty, but we want to change f to set C.A while conserving freal
    #if 1
    if ((nf.x || nf.y) && (freal >= 1-1e-6 || freal <= 1e-6)) {
        #if PRINTA
        fprintf (stderr, "|| A.S.E: f0=%0.15g alphaf=%0.15g freal=%0.15g ibm=%g\n", 
                          f0, alphaf, freal, ibm0);
        #endif

        // a. Find a value of alpha that satifies volume conservation of freal
        int itra = bisection_solver (&alpha, alphaMin, alphaMax, tcell);
        
        // b. calculate the upper value of alpha, still ensuring volume conservation
        double arangeMax = 0;
        alphaMin = alpha; alphaMax = 1;
        bisection_solver (&arangeMax, alphaMin, alphaMax, tcell);

        // c. calculate the lower value of alpha
        double arangeMin = 0;
        alphaMin = -1; alphaMax = alpha;
        bisection_solver (&arangeMin, alphaMin, alphaMax, tcell);

        // d. find the maximum value of f which satisfies fr given the alpha 
        //    range using a simplified golden-section search technique.
        golden_search (&alpha, arangeMin, alpha, arangeMax, tcell);

        return alpha;
    }
    #endif
    alphaMin = -0.75, alphaMax = 0.75;

#if 0
    double f0min = rectangle_fraction (nf, alphaMin, lhs, rhs);
    double fcalcMin = ibm0*immersed_fraction (f0min, nf, alphaMin, ns, alphas, lhs, rhs,0);
    double errorMin = fcalcMin - freal;

    f0max = rectangle_fraction (nf, alphaMax, lhs, rhs);
    fcalcMax = ibm0*immersed_fraction (f0max, nf, alphaMax, ns, alphas, lhs, rhs,0);
    errorMax = fcalcMax - freal;
#endif
    #if PRINTA
    fprintf (stderr, "|| A.S: f0=%0.15g alphaf=%0.15g freal=%0.15g ibm=%g\n", 
                      f0, alphaf, freal, ibm0);
    //fprintf(stderr, "|| f0max=%g fcalcmax=%g emax=%g f0min=%g fcalcmin=%g emin=%g\n",
    //                  f0max, fcalcMax, errorMax, f0min, fcalcMin, errorMin);
    #endif

    bisection_solver (&alpha, alphaMin, alphaMax, tcell);

    return alpha;
}


/*
this function fills fr with only the real volume fraction of f that lays 
outside of the solid immersed boundary, ibm.
*/

void real_fluid (scalar f, scalar fr)
{
    vector nf[], ns[];
    scalar alphaf[], alphas[];

    reconstruction (f, nf, alphaf);
    reconstruction (ibm, ns, alphas);

    foreach() {
        if (on_interface(ibm) && on_interface(f))
            fr[] = immersed_fraction (f[], (coord){nf.x[], nf.y[]}, alphaf[],
                                                 (coord){ns.x[], ns.y[]}, alphas[],
                                                 (coord){-0.5, -0.5, -0.5},
                                                 (coord){0.5, 0.5, 0.5}, 0);
        else
            fr[] = clamp(f[], 0., 1.);
    }
}


/* 
Some of the advected fluid gets reconstructed in the solid region of a cell.
immersed_reconstruction changes c to enforce volume/mass conservation (according
to cr) considering the immersed boundary. In other words, this function brings any
fluid otherwise reconstructed inside the solid region up in the real fluid portion
of interface cells

cr is the total real fluid in a cell after the unidimensional advection.

nf and ns are the liquid and solid normals (resp.). Likewise, alphas and alphaf
are the corresponding alpha values.
*/

void immersed_reconstruction (scalar c, const scalar cr, vector nf, scalar alphaf, 
                              vector ns, scalar alphas)
{
    foreach() {
        if (on_interface(ibm) && c[]) {

           if (cr[] == 0 || (cr[] >= 1 && c[] >= 1))
                continue;

            fprintf(stderr, "cr[] = %g\n", cr[]);

            coord nsolid = {ns.x[], ns.y[]}, nfluid = {nf.x[], nf.y[]};

//            if ((c[] <= 0 + INT_TOL || c[] >= 1 - INT_TOL) || (fabs(nfluid.x) <= 1e-10 && fabs(nfluid.y) <= 1e-10))
//                continue;
            
            double freal = immersed_fraction (c[], nfluid, alphaf[], nsolid, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);
            
            if (cr[] < freal + INT_TOL && cr[] > freal - INT_TOL) // was getting alpha convergence errors with 1e-8
                continue;

            double alpha = immersed_line_alpha (point, nfluid, alphaf[], nsolid, alphas[], cr[]);
            double c0 = c[];
            c[] = plane_volume (nfluid, alpha);
            alphaf[] = alpha;
            fprintf(stderr, "(%g, %g) c[]_before = %0.15f, c[]_after = %0.15f | cr[] = %0.15f crcalc =%0.15f\n",
                           x, y, c0, c[], cr[], freal);
        }
        //if (c[] > 1 - 1e-7) c[] = 1;
    }
#if 0
    reconstruction (c, nf, alphaf);

    scalar fr[];
    real_fluid (c, fr);

    foreach() {
        //clamp (cr[], 0., ibm[]);
        if (fr[] < cr[] - INT_TOL || fr[] > cr[] + INT_TOL) {
            double f0 = c[], fr0 = fr[];
            coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
            double freal = ibm[]*immersed_fraction (c[], nf1, alphaf[], ns1, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);

            if (freal >= cr[] - 1e-10 && freal <= cr[] + 1e-10)
                continue;

            double alpha = immersed_line_alpha (point, nf1, alphaf[], ns1, alphas[], cr[]);
            c[] = plane_volume (nf1, alpha);
            if (c[] > 1 - 1e-10) c[] = 1;

            fprintf (stderr, "clean-up (%g, %g) f0[]=%g fr0[]=%g f[]=%g fr[]=%g\n", 
                              x, y, f0, fr0, c[], cr[]);
        }
    }
    boundary ({c});
#endif
}


/*
real_volume 
*/

double real_volume (scalar f)
{
    vector nf[], ns[];
    scalar alphaf[], alphas[];

    reconstruction (f, nf, alphaf);
    reconstruction (ibm, ns, alphas);

    double volume = 0.;
    foreach(reduction(+:volume)) {
        if (on_interface(ibm) && on_interface(f))
            volume += immersed_fraction (f[], (coord){nf.x[], nf.y[]}, alphaf[],
                                              (coord){ns.x[], ns.y[]}, alphas[],
                                              (coord){-0.5, -0.5, -0.5},
                                              (coord){0.5, 0.5, 0.5}, 0) * sq(Delta);
        else
            volume += f[]*sq(Delta)*ibm[];
    }

    return volume;
}

bool is_triple_point (Point point, coord nf, coord ns);

double get_contact_angle (scalar f, scalar ibm)
{
    scalar fr_temp[];
    real_fluid (f, fr_temp);

    vector nf[], ns[];
    scalar alphaf[], alphas[];

    reconstruction (f, nf, alphaf);
    reconstruction (ibm, ns, alphas);

    double theta = 0;
    int count = 0;
    foreach(reduction(+:theta) reduction(+:count)) {
        if (on_interface(ibm) && on_interface(f) && fr_temp[]) {
            coord nf_temp = {nf.x[], nf.y[]}, ns_temp = {ns.x[], ns.y[]};
            if (is_triple_point (point, nf_temp, ns_temp)) {
                double num = nf.x[]*ns.x[] + nf.y[]*ns.y[];
                double den = distance(nf.x[], nf.y[]) * distance(ns.x[], ns.y[]);
                theta += acos (num/den);
                count++;
            }
        }
    }

    if (count > 0)
        return (theta / count)*180./pi;
    else
        return 0;
}

/*
statsf_real is a function that outputs data that you would normally get using
the statsf function but with only the portion of f that sits outside of the
solid boundary.

TODO: finish function
*/
#if 0
stats statsf_real (scalar f)
{
    vector nf0[], ns0[];
    scalar alphaf0[], alphas0[];
    reconstruction (c, nf0, alphaf0);
    reconstruction (ibm, ns0, alphas0);

    double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	        reduction(max:max) reduction(min:min)) {
        if (dv() > 0. && f[] != nodata) {
            if (on_interface(ibm) && on_interface(c)) {
                coord nft = {nf0.x[], nf0.y[]}, nst = {ns0.x[], ns0.y[]};
                coord lhs = {-0.5, -0.5}, rhs = {0.5, 0.5};
                volume += ibm[]*immersed_area(c[], nft, alphaf0[], nst, alphas0[], lhs, rhs, 0)*(sq(Delta));     
            }
            else
                volume += dv();
            sum    += dv()*f[];
            sum2   += dv()*sq(f[]);
            if (f[] > max) max = f[];
            if (f[] < min) min = f[];
        }
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;

    return s;
}
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
    ibmCells.refine = ibm_fraction_refine;
    ibmCells.prolongation = fraction_refine;

    // THESE DON'T WORK WITH AMR, EVEN IN SERIAL
    //ibmCells.refine = fraction_refine_metric;
    //ibmCells.prolongation = fraction_refine_metric;
    //ibmCells.restriction = restriction_cell_metric;

    ibm0.refine = ibm.refine = ibm_fraction_refine;
    ibm0.prolongation = ibm.prolongation = fraction_refine;

    foreach_dimension() {
        ibmFaces.x.prolongation = refine_metric_injection_x;
        ibmf.x.prolongation = ibm_face_fraction_refine_x;
        ibmFaces.x.restriction = restriction_face_metric;
    }
#endif
    restriction ({ibm, ibmf, ibmFaces, ibmCells});

    boundary(all);
}

