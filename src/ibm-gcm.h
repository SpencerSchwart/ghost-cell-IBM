#include "fractions.h"

#define BGHOSTS 2
#define IBM 1
#define LIMIT 1.e10

#undef SEPS
#define SEPS 1e-30

scalar ibm[];
scalar ibm0[];          // solid volume fraction field of previous timestep
face vector ibmf[];

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
#define on_interface(a) (a[] > 0. && a[] < 1.)
#define is_mostly_solid(a, i) (a[i] > 0. && a[i] <= 0.5)
#define is_fresh_cell(a0, a) (a0[] <= 0.5 && a[] > 0.5)


/*
These functions and macros below are used to mimic Basilisk's way of specifying
boundary condtions, e.g., u.t[embed] = dirichlet(0). For now, we only allow a
dirichlet b.c for velocity.

TODO: allow neumann b.c (navier-slip b.c?)
TODO: make other macros for p, pf, f, etc.
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
is_ghost_cell returns true if the given cell shares a face with a fluid cell,
ibm > 0.5, and the volume fraction is less than or equal to 0.5.

TODO: add check to make sure ghost cells, particularly ghost cells w/ different
      levels of refinement, dont overlap with each other.
*/

bool is_ghost_cell (Point point, scalar ibm)
{
   //return fluid_neighbor(point, ibm) && ibm[] <= 0.5 && level == depth(); // temp fix
                                                                           // for solid cells refining
   return fluid_neighbor(point, ibm) && ibm[] <= 0.5;
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


int local_to_global (Point point, coord p, double* ax, double* ay, double* az)
{
    *ax = x + p.x * Delta;
    *ay=  y + p.y * Delta;
    *az = z + p.z * Delta;

    return 0;
}


foreach_dimension()
double ibm_flux_x (Point point, scalar s, face vector mu, double * val)
{
    *val = 0.;
    if (ibm[] >= 1. || ibm[] <= 0.)
        return 0.;

    coord n, midPoint;
    double area = ibm_geometry (point, &n, &midPoint);

    double mpx, mpy, mpz;
    local_to_global(point, midPoint, &mpx, &mpy, &mpz);

    double bc = uibm_x(mpx, mpy, mpz);
    
    double coef = 0.;
    double grad = dirichlet_gradient (point, s, ibm, n, midPoint, bc, &coef);

    double mua = 0., fa = 0.;
    foreach_dimension() {
        mua += mu.x[] + mu.x[1];
        fa += fm.x[] + fm.x[1];
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
*/

/*
boundary_points is used to find the intersecting points of the fluid interface
within the area being advected.
*/

int boundary_points (coord nf, double alphaf, coord lhs, coord rhs, coord bp[2])
{
    int i = 0;
    // check on x faces first
    double dx = rhs.x - lhs.x;
    if (fabs(dx) < 1e-10)
        return 0;

    for (double xint = rhs.x; xint >= lhs.x; xint -= dx) {
        double yint = (alphaf - nf.x*xint)/nf.y;
        if (fabs(yint) <= 0.5) {
            bp[i].x = xint;
            bp[i].y = yint;
            ++i;
        }
    }

    // then check y faces
    for (double yint = -0.5; yint <= 0.5; yint += 1.) {
        double xint = (alphaf - nf.y*yint)/nf.x;
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
    pt.x = (alphas/ns.y - alphaf/nf.y) /
                  ((ns.x/ns.y) - (nf.x/nf.y));

    pt.y = (alphaf/nf.y) - (nf.x*pt.x)/nf.y;

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

double region_check (coord pc, coord nf, double alphaf, coord ns, double alphas)
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

TODO: add way to avoid a infinite loop in the while loop
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
    while (!sorted) {
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

    return fabs(area)/2;
}


/*
immersed_area calculates the area of the real fluid part of a cell.
*/

double immersed_area (coord nf, double alphaf, coord ns, double alphas, 
                  coord lhs, coord rhs)
{
    // 1. find the intersection points, pf & ps, of the fluid and solid interface
    //    with the enclosed region
    coord pf[2], ps[2];
    for (int i = 0; i < 2; ++i)
        foreach_dimension() {
            pf[i].x = nodata;
            ps[i].x = nodata;
        }
    int numpf = boundary_points(nf, alphaf, lhs, rhs, pf);
    int numps = boundary_points(ns, alphas, lhs, rhs, ps);
    
    // 2. find the intersecting point of the two interfaces (if there is one)
    coord pint = {nodata,nodata,nodata};
    int numpi = interface_intersect (nf, alphaf, ns, alphas, lhs, rhs, &pint);

    // 3. get the other points enclosing the region being advected
    coord lhst = {lhs.x,0.5}, rhsb = {rhs.x,-0.5}; 

    // 4. find which points create the polygon defining the real fluid region
    coord poly[9];
    poly[0] = pf[0], poly[1] = pf[1], poly[2] = ps[0], poly[3] = ps[1], poly[4] = pint;
    poly[5] = lhs, poly[6] = rhs, poly[7] = lhst, poly[8]= rhsb;

    int nump = 0; // # of real points
    for (int i = 0; i < 9; ++i) {
        double placement = region_check(poly[i], nf, alphaf, ns, alphas);
        if ((placement >= 0 || fabs(placement) < 1e-10) && poly[i].x != nodata) { // should be < nodata?
            nump++;
        }
        else {
            poly[i].x = nodata, poly[i].y = nodata;
        }
    }

    if (nump == 0) // we don't do anything if the region has no interface fragments
        return 0;

    coord cf[nump]; // holds real points
    int count = 0;
    for (int i = 0; i < 9; ++i) {
        if (poly[i].x != nodata && count < nump) {
            cf[count] = poly[i];
            count++;
        }
    }

    // 5. sort the real points in clockwise order
    if (lhs.x == rhs.x)
        return 0;
    sort_clockwise (nump, cf);

    // 6. use the shoelace formula to find the area
    double area = polygon_area (nump, cf);
    coord rect[4] = {lhs,rhsb,rhs,lhst};
    double areaTotal = polygon_area (4, rect);

    return area/areaTotal;

}

/*
This function calculates the fraction of a rectangle (defined by lhs and rhs)
which lies inside the interface neglecting the portion inside of the immersed boundary.

lhs and rhs are the bottom left and top right (resp.) coordinates defining the
region being advected by the split VOF advection scheme (see sweep_x in vof.h)
*/

double immersed_fraction (coord nf, double alphaf, coord ns, double alphas, 
                          coord lhs, coord rhs)
{
    double area = immersed_area(nf, alphaf, ns, alphas, lhs, rhs);
    return area;
}


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

    ibm0.refine = ibm.refine = ibm_fraction_refine;
    ibm0.refine = ibm.prolongation = fraction_refine;

    foreach_dimension() {
        ibmFaces.x.prolongation = refine_metric_injection_x;
        ibmf.x.prolongation = ibm_face_fraction_refine_x;
        ibmFaces.x.restriction = restriction_face_metric;
    }
#endif
    restriction ({ibm, ibmf, ibmFaces, ibmCells});

    boundary(all);
}

