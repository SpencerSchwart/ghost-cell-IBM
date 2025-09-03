/*
###### TWO PHASE FUNCTIONS FOR IBM ######
*/

#define PRINTA 0
#define BI_TOL 1e-9 // TODO: standardize tolerance for root solvers

#include "fractions.h"

(const) scalar contact_angle;

extern scalar f;

typedef struct proj_coord {
    coord og;
    double x, y;
} proj_coord;

typedef struct triangle {
    coord p[3];
} triangle;

typedef struct plane {
    coord n;
    double alpha;

    coord * p;      // stores all of the points along the plane
    int psize;      // size of p array / number of points

    triangle * tri; // stores coordinates for each triangle
    int tsize;      // # of triangles
} plane;

double get_percent_error (double a0, double a1)
{
    return (a1 - a0)/(a0+SEPS) * 100;
}

int percent_error_tol (double a0, double a1, double TOL = VTOL)
{
    return fabs(get_percent_error(a0, a1)) <= TOL;
}

coord get_centroid (int size, coord * np)
{
    coord centroid = {0,0,0};
    for (int i = 0; i < size; ++i) {
        foreach_dimension()
            centroid.x += np[i].x;
    }
    
    foreach_dimension()
        centroid.x /= size + SEPS;

    return centroid;
}

int near_id_cell (Point point, scalar id)
{
    foreach_neighbor(1) { // only 3x3x3 stencil checked
        if (id[] > 0)
            return 1;
    }
    return 0;
}

trace
void reconstruction_real (const scalar c, vector n, scalar alpha, scalar id)
{
    foreach() {
        if (c[] <= 0. || c[] >= 1.) {
            alpha[] = 0.;
            foreach_dimension()
                n.x[] = 0.;
        }
        else if (near_id_cell(point, id)) {
            coord m = interface_normal (point, c);
            foreach_dimension()
                n.x[] = m.x;
            alpha[] = plane_alpha(c[], m);
        }
    }

#if TREE // is this redundant?
    foreach_dimension()
        n.x.refine = n.x.prolongation = refine_injection;

    alpha.n = n;
    alpha.refine = alpha.prolongation = alpha_refine;
#endif
}

void print_plane_points(int pcount, plane * planes)
{
    for (int i = 0; i < pcount; ++i) {
        fprintf(stderr, "\nplane #%d n={%g, %g, %g} alpha=%g\n", 
        i+1, planes[i].n.x, planes[i].n.y, planes[i].n.z, planes[i].alpha);

        for (int j = 0; j < planes[i].psize; ++j) {
            fprintf(stderr, "point #%d p = {%g, %g, %g}\n",
                j+1, planes[i].p[j].x, planes[i].p[j].y, planes[i].p[j].z);
        }
    }
}

void print_plane_triangles(int pcount, plane * planes)
{
    for (int i = 0; i < pcount; ++i) {
        fprintf(stderr, "\nplane #%d n={%g, %g, %g} alpha=%g\n", 
        i+1, planes[i].n.x, planes[i].n.y, planes[i].n.z, planes[i].alpha);

        for (int j = 0; j < planes[i].tsize; ++j) {

            fprintf (stderr, "triangle #%d\n", j+1);
            for (int k = 0; k < 3; ++k) {
                fprintf(stderr, "point #%d p = {%g, %g, %g}\n",
                    k+1, planes[i].tri[j].p[k].x, planes[i].tri[j].p[k].y, planes[i].tri[j].p[k].z);
            }
        }
    }
}

void print_coord_list (int size, coord list[size])
{
    for (int i = 0; i < size; ++i) {
        fprintf(stderr, "point #%d = {%0.15g, %0.15g, %0.15g}\n", i+1, list[i].x, list[i].y, list[i].z);
    }
}


// returns p given t for the parametric equation
coord plane_p0 (plane pl1, plane pl2, coord d)
{
    coord n1 = pl1.n, n2 = pl2.n;
    double alpha1 = pl1.alpha, alpha2 = pl2.alpha;

    int drop = (fabs(d.x) >= fabs(d.y) && fabs(d.x) >= fabs(d.z)) ? 0 :
               (fabs(d.y) >= fabs(d.z)) ? 1 : 2;

    double a1, a2, b1, b2;

    switch (drop) {
        case 0: // x = t = 0
            a1 = n1.y; b1 = n1.z; a2 = n2.y; b2 = n2.z; break;
        case 1: // y = 0
            a1 = n1.x; b1 = n1.z; a2 = n2.x; b2 = n2.z; break;
        case 2: // z = 0
            a1 = n1.x; b1 = n1.y; a2 = n2.x; b2 = n2.y; break;
    }

    double det = a1*b2 - a2*b1;

    //assert (fabs(det) > 1e-14);

    double u = (alpha1 * b2 - alpha2 * b1) / (det + SEPS);
    double v = (alpha2 * a1 - alpha1 * a2) / (det + SEPS);

    coord p0;
    switch (drop) {
        case 0:
            p0.x = 0; p0.y = u; p0.z = v; break;
        case 1:
            p0.x = u; p0.y = 0; p0.z = v; break;
        case 2:
            p0.x = u; p0.y = v; p0.z = 0; break;
    }

    return p0;
}

double plane_t (coord p0, coord d, coord face)
{
    int type = face.x? 0: face.y? 1: 2;

    switch(type) {
        case 0:
            if (fabs(d.x) < 1e-14) return HUGE;
            return (face.x - p0.x) / d.x;
        case 1:
            if (fabs(d.y) < 1e-14) return HUGE;
                return (face.y - p0.y) / d.y;
        case 2:
            if (fabs(d.z) < 1e-14) return HUGE;
                return (face.z - p0.z) / d.z;                
    }
    return HUGE;
}

int bounds_check (coord p, double ufdx = 0.5)
{
    return p.x <= ufdx+VTOL && fabs(p.x) <= 0.5+VTOL && fabs(p.y) <= 0.5+VTOL && fabs(p.z) <= 0.5+VTOL;
}

int point_is_unique (int size, coord pint[size], coord p)
{
    for (int i = 0; i < size; ++i)
        if (approx_equal(pint[i], p, VTOL))
            return 0;
    return 1;
}

int plane_is_unique (int size, plane pint[size], plane p)
{
    for (int i = 0; i < size; ++i)
        if (approx_equal(pint[i].n, p.n, VTOL) && 
            approx_equal_double(pint[i].alpha, p.alpha, VTOL))
            return 0;
    return 1;
}

const coord cubeface[6] = {
        { 0.5, 0.0, 0.0},
        {-0.5, 0.0, 0.0},
        { 0.0, 0.5, 0.0},
        { 0.0,-0.5, 0.0},
        { 0.0, 0.0, 0.5},
        { 0.0, 0.0,-0.5}
};

// fills pint with the intersection points of the provided two planes if
// they exist and are located within the bounded region of the cell
int plane_intersect (plane pl1, plane pl2, plane padv, coord pint[2], int print = 0)
{
    if (determinant(pl1.n, pl2.n) <= 1e-15) {
        #if PRINTA
        fprintf(stderr, "planes are parallel\n");
        #endif
        return 0;   // planes are parallel
    }

    coord d = cross_product(pl1.n, pl2.n);
    coord p0 = plane_p0 (pl1, pl2, d);

#if PRINTA
if (print) {
    fprintf(stderr, "d = \n");
    print_coord_list(1, &d);

    fprintf(stderr, "p0 = \n");
    print_coord_list(1, &p0);
}
#endif

    coord cubeface_copy[6];
    memcpy(cubeface_copy, cubeface, sizeof(cubeface));
    cubeface_copy[0].x = padv.alpha; // change +x plane to the advection plane

    int count = 0;
    for (int i = 0; i < 6; ++i) {
        double t = plane_t (p0, d, cubeface_copy[i]);

        if (t == HUGE)
            continue;

        coord pit = {p0.x + t*d.x, p0.y + t*d.y, p0.z + t*d.z};

#if PRINTA
if (print) {
        fprintf(stderr, "pi = \n");
        print_coord_list(1, &pit);
}
#endif

        //if (approx_equal (pit, p0, VTOL))
        //    continue;

#if PRINTA
if (print) {
        fprintf(stderr, "point != p0\n");

        fprintf(stderr, "bounds_check = %d, unique = %d\n", 
            bounds_check(pit, padv.alpha), point_is_unique(2, pint, pit));
}
#endif

        if (bounds_check(pit, padv.alpha) && point_is_unique(2, pint, pit)) {
            //assert(count < 2);
            if (count > 2)
                fprintf(stderr, "WARNING: count > 2 in plane_intersect!\n");
            else {
                pint[count] = pit;
                count++;
            }
        }
    }

    return count;
}

static coord cube_vertices[8] = {
        { 0.5, 0.5, 0.5},
        { 0.5,-0.5, 0.5},
        { 0.5, 0.5,-0.5},
        { 0.5,-0.5,-0.5},
        {-0.5, 0.5, 0.5},
        {-0.5, 0.5,-0.5},
        {-0.5,-0.5, 0.5},
        {-0.5,-0.5,-0.5}
};

/*
This function returns true if there exists a triple points within the cell, 
i.e. the liquid interface intersects the solid one, given the normals of the
liquid and solid interfaces.
*/

bool is_triple_point (Point point, coord nf, coord ns)
{
    if (!(on_interface(ibm)) || !(on_interface(f)))
        return false;
    if ((ns.x == 0 && ns.y == 0 && ns.z == 0) || (nf.x == 0 && nf.y == 0 && nf.z == 0))
        return false;

    double alphas = plane_alpha (ibm[], ns);

    double alphaf = plane_alpha (f[], nf);

    #if 0
    double intercept = ((alphas/(ns.y+SEPS)) - (alphaf/(nf.y+SEPS))) /
                       ((ns.x/(ns.y+SEPS)) - (nf.x/(nf.y+SEPS)) + SEPS);
    #endif

#if dimension == 2
    coord pt_int = {0,0};
    pt_int.x = ((alphas/(ns.y+SEPS)) - (alphaf/(nf.y+SEPS))) /
                  ((ns.x/(ns.y+SEPS)) - (nf.x/(nf.y+SEPS)) + SEPS);
    pt_int.y = (alphaf/(nf.y + SEPS)) - (nf.x*pt_int.x)/(nf.y + SEPS);

    return (fabs(pt_int.x) <= 0.5 && fabs(pt_int.y) <= 0.5);
#else
    plane plf = {nf, alphaf};
    plane pls = {ns, alphas};
    plane padv = {{1,0,0}, 0.5};
    coord pint[2];

    return plane_intersect (plf, pls, padv, pint);
#endif
}


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

TODO: do we need another case when ns.y || nf.y = 0? or does this code
      handle everything ok?
*/

int interface_intersect (coord nf, double alphaf, coord ns, double alphas,
                         coord lhs, coord rhs, coord * pint = NULL)
{
#if dimension == 2
    coord pt;
    pt.x = (alphas/(ns.y + SEPS) - alphaf/(nf.y + SEPS)) /
                  ((ns.x/(ns.y + SEPS)) - (nf.x/(nf.y + SEPS)) + SEPS);

    pt.y = (alphaf/(nf.y + SEPS)) - (nf.x*pt.x)/(nf.y + SEPS);

    foreach_dimension() {
        if (pt.x < lhs.x || pt.x > rhs.x) {
            return 0;
        }
    }
    if (pint)
        *pint = pt;
#endif
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

double region_check2(plane plf, plane pls, coord pc)
{
#if dimension == 2
    double fluid = plf.alpha - plf.n.x*pc.x - plf.n.y*pc.y;
    double solid = pls.alpha - pls.n.x*pc.x - pls.n.y*pc.y;
#else
    double fluid = plf.alpha - plf.n.x*pc.x - plf.n.y*pc.y - plf.n.z*pc.z;
    double solid = pls.alpha - pls.n.x*pc.x - pls.n.y*pc.y - pls.n.z*pc.z;
#endif
    //return difference(fluid,-solid); //- because solid has inward pointing normal
    //return intersection(fluid,solid); // these options are equivalent
    return min(fluid,solid);
}

// + = inside, - = outside, 0 = on interface
double region_check (plane pl, coord pc)
{
#if dimension == 2
    return pl.alpha - pl.n.x*pc.x - pl.n.y*pc.y;
#else
    return pl.alpha - pl.n.x*pc.x - pl.n.y*pc.y  - pl.n.z*pc.z;
#endif
}

int vertices_region(plane pls, plane plf, plane padv, coord pv[8], int print = 0)
{
    coord cube_vertices_copy[8];
    memcpy(cube_vertices_copy, cube_vertices, sizeof(cube_vertices));

    for (int i = 0; i < 8; ++i) {
        if (approx_equal_double (cube_vertices_copy[i].x, 0.5))
            cube_vertices_copy[i].x = padv.alpha;
    }

    int count = 0;
    for (int i = 0; i < 8; ++i) {
        double placement = region_check2(plf, pls, cube_vertices_copy[i]);
        if ((placement >= 0 || fabs(placement) < 1e-6) && bounds_check(cube_vertices_copy[i], padv.alpha))
            pv[i] = cube_vertices_copy[i], count++;
        else
            foreach_dimension()
                pv[i].x = HUGE;
    }
    return count;
}

int fill_faces (plane plf, plane pls, plane padv, int nump, const coord tp[nump], plane** planes)
{
    //plane cellpx = {{ 1, 0, 0}, 0.5};
    //plane cellpx = padv;
    plane cellnx = {{-1, 0, 0}, 0.5};
    plane cellpy = {{ 0, 1, 0}, 0.5};
    plane cellny = {{ 0,-1, 0}, 0.5};
    plane cellpz = {{ 0, 0, 1}, 0.5};
    plane cellnz = {{ 0, 0,-1}, 0.5};

    int numplanes = 6;

    plane allPlanes[8] = {
        plf, cellnx, cellpy, cellny, cellpz, cellnz
    };

    // edge case, for when advection plane coincides with another plane
    if (plane_is_unique(numplanes, allPlanes, padv)) {
        allPlanes[numplanes] = padv;
        ++numplanes;
    }

    // edge case, for when solid plane coincides with another plane (typically liquid plane)
    if (plane_is_unique(numplanes, allPlanes, pls)) {
        allPlanes[numplanes] = pls;
        ++numplanes;
    }

    int planeCount = 0;
    for (int i = 0; i < numplanes; ++i) { // cycle through each plane;

        int pcount = 0;
        for (int j = 0; j < nump; ++j) { // check each point
            double placement = region_check(allPlanes[i], tp[j]);
            if (fabs(placement) <= 1e-10) // point is on face
                pcount++;
        }

        if (pcount < 3) {
            allPlanes[i].p = NULL;
            allPlanes[i].psize = 0;
            continue;
        }
        
        assert (pcount >= 3);
        
        allPlanes[i].p = malloc((pcount * sizeof(coord)));
        allPlanes[i].psize = pcount;

        for (int j = 0, jtrue = 0; j < nump; ++j) {
            double placement = region_check(allPlanes[i], tp[j]);
            if (fabs(placement) <= 1e-10) {
                allPlanes[i].p[jtrue] = tp[j];
                jtrue++;
            }
        }

        planeCount++;
    }

    *planes = malloc(planeCount * sizeof(plane));

    for (int i = 0, itrue = 0; i < numplanes; ++i) {
        if (allPlanes[i].psize > 0) {
            (*planes)[itrue] = allPlanes[i];
            ++itrue;
        }
    }

    return planeCount;
}

/*
The next few functions are to allow portablity for qsort_r
*/

typedef int (*qsort_cmp_r)(const void *a, const void *b, void *arg);

static void *qsort_r_arg;

static int qsort_r_trampoline(const void *a, const void *b) {
    extern qsort_cmp_r qsort_r_func;
    return qsort_r_func(a, b, qsort_r_arg);
}

qsort_cmp_r qsort_r_func;

void qsort_r_fallback(void *base, size_t nmemb, size_t size,
                      qsort_cmp_r compar, void *arg) {
    qsort_r_func = compar;
    qsort_r_arg = arg;
    qsort(base, nmemb, size, qsort_r_trampoline);
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

int is_behind (const void *pa, const void *pb, void *center)
{
    const coord *a = pa, *b = pb, *c = center;
    double atheta = atan2 (a->y - c->y, a->x - c->x);
    double btheta = atan2 (b->y - c->y, b->x - c->x);

    if (atheta < 0) atheta += 2*pi;
    if (btheta < 0) btheta += 2*pi;

    if (atheta > btheta) return -1;
    if (atheta < btheta) return 1;

    // a and b lay along the same line, so det = 0. Instead, we check to see
    // if a is further from the center than b.
    double dista = sq(a->x - c->x) + sq(a->y - c->y);
    double distb = sq(b->x - c->x) + sq(b->y - c->y);

    return (dista > distb) ? -1: (dista < distb);
}

int compare_projected(const void *a, const void *b) {
    const proj_coord *pa = a;
    const proj_coord *pb = b;
    double angle_a = atan2(pa->y, pa->x);
    double angle_b = atan2(pb->y, pb->x);
    return (angle_a < angle_b) ? -1 : (angle_a > angle_b);
}

coord face_normal(coord a, coord b, coord c)
{
    coord ab = {b.x - a.x, b.y - a.y, b.z - a.z};
    coord ac = {c.x - a.x, c.y - a.y, c.z - a.z};
    coord n = cross_product (ab, ac);

    double norm = 0;
    foreach_dimension()
        norm += n.x;

    foreach_dimension()
        n.x /= norm + SEPS;
    
    return n;  
}


/*
sort_clockwise sorts a list of coordinates, provided in cf w/nump points, in
clockwise order (or counter-clockwise if y-advection).

TODO: make custome qsort_r function to allow compatiblity with Mac compilers
*/

void sort_clockwise (int nump, coord cf[nump], int print = 0)
{
    double xsum = 0, ysum = 0;
    for (int i = 0; i < nump; ++i) {
        xsum += cf[i].x;
        ysum += cf[i].y;
    }
    coord pc = {xsum/nump, ysum/nump}; // center coordinate is average of all points
   
    if (print)
        fprintf(stderr, "SORTING: pc={%g, %g} nump=%d\n", pc.x, pc.y, nump);
    qsort_r_fallback (cf, nump, sizeof(coord), is_behind, &pc);
}

void sort_clockwise3 (int nump, coord cf[nump], coord n = {0,0,0})
{
    assert (nump >= 3);

    coord a = cf[0], b = cf[1];
    if (!n.x && !n.y && !n.z) {
        coord c = cf[2];
        n = face_normal(a, b, c); // if n isn't provided
    }
    normalize2(&n);

    coord u = (coord){b.x - a.x, b.y - a.y, b.z - a.z};
    coord v = cross_product(n, u);

    normalize2(&u); normalize2(&v);

    coord centroid = get_centroid(nump, cf);

    proj_coord cfp[nump];

    for (int i = 0; i < nump; ++i) {
        double px = (cf[i].x - centroid.x)*u.x + 
                    (cf[i].y - centroid.y)*u.y +
                    (cf[i].z - centroid.z)*u.z;
        double py = (cf[i].x - centroid.x)*v.x + 
                    (cf[i].y - centroid.y)*v.y +
                    (cf[i].z - centroid.z)*v.z;
        cfp[i].og = cf[i];
        cfp[i].x = px, cfp[i].y = py;
    }

    qsort(cfp, nump, sizeof(proj_coord), compare_projected);

    for (int i = 0; i < nump; ++i) {
        cf[i] = cfp[i].og;

        #if PRINTA
        //fprintf(stderr, "cf = {%g, %g, %g} cfp = {%g, %g}\n",
        //cf[i].x, cf[i].y, cf[i].z, cfp[i].x, cfp[i].y);
        #endif
    }
    #if PRINTA
    //fprintf(stderr, "\n");
    #endif
}

void triangulate_planes(int plcount, plane * planes)
{
    for (int i = 0; i < plcount; ++i) {

        assert (planes[i].psize >= 3);

        planes[i].tri = malloc ((planes[i].psize - 2) * sizeof(triangle));
        planes[i].tsize = planes[i].psize - 2;

        for (int j = 1, tricount = 0; j < planes[i].psize - 1; ++j) {
            assert (tricount < (planes)[i].tsize);
            planes[i].tri[tricount] = (triangle){{planes[i].p[0], planes[i].p[j], planes[i].p[j+1]}};
            tricount++;
        }
    }
}

double tetrahedron_volume(triangle tri, coord pc)
{
    coord a = {tri.p[0].x - pc.x, tri.p[0].y - pc.y, tri.p[0].z - pc.z};
    coord b = {tri.p[1].x - pc.x, tri.p[1].y - pc.y, tri.p[1].z - pc.z};
    coord c = {tri.p[2].x - pc.x, tri.p[2].y - pc.y, tri.p[2].z - pc.z};

    return 1./6. * fabs((a.x*(b.y*c.z - c.y*b.z) +
                         a.y*(b.z*c.x - c.z*b.x) +
                         a.z*(b.x*c.y - c.x*b.y)));
}

double polyhedron_volume(int pcount, plane * planes, coord pc)
{
    double volume = 0;
    for (int i = 0; i < pcount; ++i)
        for (int j = 0; j < planes[i].tsize; ++j) {
            volume += tetrahedron_volume(planes[i].tri[j], pc);
        }
    
    return volume;
} 

double rectangular_prism_volume (coord bc, coord tc)
{
    double l = tc.x - bc.x;
    double w = tc.y - bc.y;
    double h = tc.z - bc.z;

    return l*w*h;
}

int make_list_unique(coord** tp, int tsize, int sizes[tsize], coord* set[tsize], 
                     plane plf, plane pls)
{
    coord unique[34];
    int ucount = 0;

    for (int i = 0; i < tsize; ++i) {
        for (int j = 0; j < sizes[i]; ++j) {
            double placement = region_check2(plf, pls, set[i][j]);
            if ((placement > 0 || fabs(placement) <= 1e-6) && set[i][j].x != HUGE && 
                point_is_unique(ucount, unique, set[i][j])) {
                unique[ucount] = set[i][j];
                ucount++;
            }
        }
    }

    *tp = malloc (ucount * sizeof(coord));

    int ucount1 = 0;
    for (int i = 0; i < tsize; ++i) {
        for (int j = 0; j < sizes[i]; ++j) {
            double placement = region_check2(plf, pls, set[i][j]);
            if ((placement > 0 || fabs(placement) <= 1e-6) && set[i][j].x != HUGE && 
                point_is_unique(ucount1, *tp, set[i][j])) {
                (*tp)[ucount1] = set[i][j];
                ucount1++;
            }
        }
    }

    return ucount;
}                         

int remove_invalid_points(int size, coord ps[size], plane fluid, plane solid, plane padv)
{
    int fcount = 0;
    for (int i = 0; i < size; ++i) {
        double placement = region_check2(fluid, solid, ps[i]);
        if ((placement >= 0 || fabs(placement) < 1e-6) && bounds_check(ps[i], padv.alpha))
            fcount++;
        else
            ps[i].x = nodata, ps[i].y = nodata, ps[i].z = nodata;
    }
    return fcount;
}

double fit_volume (double advVolume, coord ns, double alphas, double ufdt);

#if dimension == 3
double immersed_volume (double c, plane plf, plane pls, coord lhs, coord rhs, 
                        double advVolume = 0, int print = 0)
{
    if (lhs.x == rhs.x || c <= 0)
        return 0;

    if (rhs.x < 0.5 && advVolume) {
        rhs.x = fit_volume (advVolume, pls.n, pls.alpha, rhs.x);
    }

    double ufdt = rhs.x; // advection plane is always the +x plane
    const plane padv = {{1,0,0}, ufdt};

    coord pint[2] = {{nodata, nodata, nodata}, {nodata, nodata, nodata}};
    plane_intersect(plf, pls, padv, pint, print);

    if (print)
        print_coord_list (2, pint);

    coord fp[12];
    int fcount0 = facets(plf.n, plf.alpha, fp, 1);
    remove_invalid_points(fcount0, fp, plf, pls, padv);

    coord sp[12];
    int scount0 = facets(pls.n, pls.alpha, sp, 1);

    remove_invalid_points(scount0, sp, plf, pls, padv);

    coord vp[8];
    vertices_region(pls, plf, padv, vp, print);

    if (print)
        print_coord_list (8, vp);

    coord pfint[2] = {{nodata, nodata, nodata}, {nodata, nodata, nodata}};
    plane_intersect(plf, padv, padv, pfint, print);

    if (print)
        print_coord_list (2, pfint);

    coord psint[2] = {{nodata, nodata, nodata}, {nodata, nodata, nodata}};
    plane_intersect(pls, padv, padv, psint, print);

    if (print)
        print_coord_list (2, psint);

    coord* totalSet[6] = {pint, fp, sp, vp, pfint, psint}, *tp = NULL;
    int totalSetSize[6] = {2, fcount0, scount0, 8, 2, 2};
    int rcount = make_list_unique(&tp, 6, totalSetSize, totalSet, plf, pls);

    plane * planes;
    int planeCount = fill_faces(plf, pls, padv, rcount, tp, &planes);

    for (int i = 0; i < planeCount; ++i)
        sort_clockwise3(planes[i].psize, planes[i].p, planes[i].n);

    triangulate_planes(planeCount, planes);
    coord centroid = get_centroid (rcount, tp);

    double realVolume = polyhedron_volume (planeCount, planes, centroid);
    double totalVolume = rectangular_prism_volume (lhs, rhs); 
    double liquidVolume = rectangle_fraction (pls.n, pls.alpha, lhs, rhs);

    // The advected volume from the argument may be larger than what is possible 
    // (small cell problem). Therefore, we only calculate it again if it is not provided.
    if (!advVolume)
        advVolume = liquidVolume*totalVolume;

    double vf = clamp (realVolume/advVolume, 0., 1.);

    #if PRINTA
    if (advVolume && !approx_equal_double(liquidVolume*totalVolume, advVolume))
        fprintf(stderr, "WARNING: advected volumes do not equal eachother! given = %g calc = %g\n",
            advVolume, liquidVolume*totalVolume);
    #endif

//#if PRINTA
    if (print) {
      fprintf(stderr, "fluid = ({%0.15g, %0.15g, %0.15g}, %0.15g)\nsolid = ({%g, %g, %g}, %g) adv=%g\n\n",
            plf.n.x, plf.n.y, plf.n.z, plf.alpha, pls.n.x, pls.n.y, pls.n.z, pls.alpha, ufdt);
        print_coord_list(rcount, tp);
        print_plane_points(planeCount, planes);
        //print_plane_triangles(planeCount, planes);
        fprintf(stderr, "real volume = %g, total volume = %g, liquid volume = %g, advVolume=%g, vf=%g\n",
            realVolume, totalVolume, liquidVolume, advVolume, vf);
    }
//#endif

    for (int i = 0; i < planeCount; ++i) {
        if (planes[i].p) free (planes[i].p);
        if (planes[i].tri) free (planes[i].tri);
    }
    free (planes);
    free (tp);

    return vf;
}
#endif


/*
polygon_area calculates the area enclosed by a list of points (in cw or ccw order)
using the shoelace formula.
*/

double polygon_area (int nump, coord cf[nump])
{
    double area = 0;
    for (int i = 0; i < nump; ++i) {
        int next = i + 1 < nump? i + 1: 0; // to close the shape
        area += cf[i].x*cf[next].y - cf[next].x*cf[i].y;
    }

    return fabs(area)/2.;
}

double rectangle_area (coord bp, coord tp)
{
    return (tp.x - bp.x) * (tp.y - bp.y); // width * height
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


/*
immersed_area calculates the area of the real fluid part of a cell given a 
bounding box (defined in lhs and rhs)
*/

double immersed_area (double c, coord nf, double alphaf, coord ns, double alphas, 
                      coord lhs, coord rhs, double advVolume = 0, int print = 0)
{
    if (lhs.x == rhs.x)
        return 0;

    // 1. calculate the area of the real region within the advected region
    coord lhst = {lhs.x,0.5}, rhsb = {rhs.x,-0.5}; 
    coord rect[4] = {lhs,rhsb,rhs,lhst};
    double areaTotal = polygon_area (4, rect); // total area being considered for advection (uf*dt*h)
    double areaLiquid = rectangle_fraction (ns, alphas, lhs, rhs); // volume fraction that isn't solid

    // TODO: temp fix when flux area > fluid area (small cell problem)
    #if 0
    if (alphas - ns.x*rhs.x - ns.y*rhs.y < 0 && alphas - ns.x*rhsb.x - ns.y*rhsb.y < 0)
        areaLiquid = 1.;
    #endif

    //#if PRINTA
    double cvy = line_intersect (alphas, ns, x = -0.5); // intercept of solid interface on left face
    //#endif

    #if 1 // temporarily off
    double rhsx = 0;
    if (rhs.x < 0.5 && areaLiquid < 1 && advVolume) {
        cvy = clamp(cvy, -0.5, 0.5);
        rhsx = fit_volume (advVolume, ns, alphas, rhs.x);
        rhs.x = rhsx;
        rhsb.x = rhsx;
        coord rect0[4] = {lhs,rhsb,rhs,lhst};
        areaTotal = polygon_area (4, rect0);
        areaLiquid = rectangle_fraction (ns, alphas, lhs, rhs);
    }
    #endif

    double areaAdv = 0;

    #if 1
    if (!advVolume)
        areaAdv = areaTotal*areaLiquid; // = uf*dt*ibmf*dx
    else
        areaAdv = advVolume;
    #else
        areaAdv = areaTotal*areaLiquid; // = uf*dt*ibmf*dx
    #endif

    #if PRINTA
    if (advVolume && !approx_equal_double(areaTotal*areaLiquid, advVolume))
        fprintf(stderr, "WARNING: advected volumes do not equal eachother! given = %g calc = %g\n",
            advVolume, areaTotal*areaLiquid);
    #endif

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
        double placement = region_check2((plane){nf, alphaf}, (plane){ns, alphas}, poly[i]);
        if ((placement >= 0 || fabs(placement) < VTOL) && poly[i].x != nodata) // should be < nodata?
            nump++;
        else 
            poly[i].x = nodata, poly[i].y = nodata;
    }

    if (print == 1) {
        fprintf(stderr, "||  pf1=(%g,%g) pf2(%g,%g) ps1=(%g,%g) ps2=(%g,%g) pint=(%g,%g)\n",
                        pf[0].x, pf[0].y, pf[1].x, pf[1].y, ps[0].x, ps[0].y, ps[1].x, 
                        ps[1].y, pint.x, pint.y);
   
        fprintf(stderr, "|| p0=(%g,%g) p1=(%g,%g) p2=(%g,%g) p3=(%g,%g) p4=(%g,%g)"
                        "   p5=(%g,%g) p6=(%g,%g) p7=(%g,%g) p8=(%g,%g)\n", 
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
    sort_clockwise (nump, cf, print);

    // 6. use the shoelace formula to find the area
    double area = polygon_area (nump, cf);

    #if PRINTA
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
    #endif

    double vf = clamp(area/(areaAdv+SEPS), 0., 1.);
    return vf;
}

/*
This function calculates the fraction of a rectangle (defined by lhs and rhs)
which lies inside the liquid interface neglecting the portion inside of the 
immersed boundary.

lhs and rhs are the bottom left and top right (resp.) coordinates defining the
region being advected by the split VOF advection scheme (see sweep_x in vof.h)
*/

trace
double immersed_fraction (double c, coord nf, double alphaf, coord ns, double alphas, 
                          coord lhs, coord rhs, double advVolume = 0, int print = 0)
{
    if (c <= 0)
        return 0;

#if dimension == 2
    return immersed_area(c, nf, alphaf, ns, alphas, lhs, rhs, advVolume, print);
#else
    plane liquid = {nf, alphaf};
    plane solid = {ns, alphas};

    if (print)
    fprintf(stderr, "Immersed_fraction: nf={%0.15g, %0.15g, %0.15g} alphaf=%0.15g ns={%g, %g, %g} alphas=%g\n",
                nf.x, nf.y, nf.z, alphaf, ns.x, ns.y, ns.z, alphas);


    return immersed_volume(c, liquid, solid, lhs, rhs, advVolume, print);
#endif
}



/*
immersed_line_alpha calculates the alpha value that conserves the volume of
real fluid, given in freal. We find the root of a function using the iterative
bisection method to obtain alpha.

The solver converges after about 20 iterations if tolerance = 1e-7
*/


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

// generic root solver using the bisection method
int rsolver_bisection (double* a, double amin, double amax, const void* data, double (*func)(const void*, double),
                       double tolerance = BI_TOL, int maxitr = 50, int print = 0)
{
    double errmax = (*func)(data, amax);
    if (fabs(errmax) < tolerance) { // problematic for receding, hydrophobic CAs?
        *a = amax;
        return 0;
    }

    double b = 0, error = HUGE;
    int itr = 0;
    while (fabs(error) > tolerance && itr < maxitr) {
        b = (amin + amax)/2.;   // bisection
        error = (*func)(data, b);
        if (sign2(error) == sign2(errmax)) {
            amax = b;
            errmax = error;
        }
        else
            amin = b;
        itr++;
    }
        
    if (itr == maxitr) {
        fprintf(stderr, "WARNING: alpha  solver does not converge after"
                        " maximum iteration (%d), error = %g\n", maxitr, error);
    }
    *a = b;
    return itr;
}

// generic root solver using Brent's method
int rsolver_brent (double* result, double a, double b, const void* data, double (*func)(const void*, double),
                   double tolerance = BI_TOL, int maxitr = 50, int print = 0, int warning = 1)
{
    double fa = (*func)(data, a);
    if (fabs(fa) < tolerance) {
        *result = a;
        return 0;
    }

    double fb = (*func)(data, b);
    if (fabs(fb) < tolerance) {
        *result = b;
        return 0;
    }

    if (fa * fb >= 0) {
        if (warning)
            fprintf(stderr, "WARNING: range in Brent's solver does not contain root!\n");

    #if PRINTA
        fprintf(stderr, "    f(%g) = %0.15g, f(%g) = %0.15g\n",
            a, fa, b, fb);
    #endif

        return -1;
    }

    if (fabs(fa) < fabs(fb)) {
        swap(double, a, b);
        swap(double, fa, fb);
    }

    double c = a, fc = fa, d = b - a; 
    int iter = 0, mflag = 1;
    //while (fabs(fb) > tolerance && fabs(fa) > tolerance && iter < maxitr) {
    while (fabs(b - a) > tolerance && fabs(fb) > tolerance && fabs(fa) > tolerance && iter < maxitr) {

        double s;
        if (!approx_equal_double(fa, fc) && !approx_equal_double(fb, fc)) {
            s = (a*fb*fc)/((fa-fb)*(fa-fc)) + 
                (b*fa*fc)/((fb-fa)*(fb-fc)) + 
                (c*fa*fb)/((fc-fa)*(fc-fb));
        }
        else // secant
            s = b - fb*(b-a)/(fb-fa + SEPS);

        int cond1 = s <= fmin(a, b) || s >= fmax(a, b);
        //int cond1 = (s < (3*a + b)/4.) || (s > b); // assuming a < b and |fa| > |fb|
        int cond2 =  mflag && fabs(s - b) >= fabs(b-c)*0.5;
        int cond3 = !mflag && fabs(s - b) >= fabs(c-d)*0.5;
        int cond4 =  mflag && fabs(b - c) < tolerance;
        int cond5 = !mflag && fabs(c - d) < tolerance;
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = (a + b) * 0.5;
            mflag = 1;
        }
        else
            mflag = 0;

        double fs = (*func)(data, s);
        d = c;        
        c = b; fc = fb;

        if (fa*fs < 0)
            b = s, fb = fs;
        else
            a = s, fa = fs;

        if (fabs(fa) < fabs(fb)) {
            swap(double, a, b);
            swap(double, fa, fb);
        }
        if (fabs(fb) < tolerance || fabs(fa) < tolerance)
            break;

        ++iter;
    }

    if (iter == maxitr && warning) 
        fprintf(stderr, "WARNING: alpha  solver does not converge after"
                        " maximum iteration (%d), error = %g\n", maxitr, fb);

    *result = (fabs(fa) < fabs(fb)) ? a : b;
    return iter;
}

/** 
adjust x coordinate (rhsx) to conserve area
    lhsb = left-hand-side bottom, rhst = right-hand-side top
*/

double fit_volume_error (const void* data, double newx)
{
    const tripoint tcell = (*(tripoint*)data);
    
    static const coord lhs = {-0.5,-0.5,-0.5};
    coord rhs = {newx, 0.5, 0.5};

#if dimension == 2
    double vtotal = rectangle_area (lhs, rhs);
#else
    double vtotal = rectangular_prism_volume (lhs, rhs);
#endif
    double vfrac = rectangle_fraction (tcell.ns, tcell.alphas, lhs, rhs);
    double vreal = vtotal * vfrac;

#if PRINTA && 0
    fprintf(stderr, "newx=%0.15g vtotal=%0.15g vfrac=%0.15g vreal=%0.15g tcell.fr=%0.15g err=%g\n",
        newx, vtotal, vfrac, vreal, tcell.fr, (vreal - tcell.fr));
#endif

    // TODO: divide by dv() or something to allow for scalability
    return (vreal - tcell.fr); // in this case, tcell.fr holds advVolume
}

// adjust ufdt to match the advected volume with the real advection volume
double fit_volume (double advVolume, coord ns, double alphas, double ufdt)
{
    if (advVolume < 1e-9)
        return ufdt;

    double oldx = ufdt;
    
    tripoint tcell;
    tcell.ns = ns; tcell.alphas = alphas; tcell.fr = advVolume;

    double xmin = -0.5, xmax = 0.5, newx = HUGE;

    int maxitr = 30;
    int itr = rsolver_brent (&newx, xmin, xmax, &tcell, fit_volume_error, 
                             tolerance = 1e-14, maxitr = maxitr, warning = 0);
    
    #if 1
    bool smallcell = false;
    if (itr == -1) { // small cell problem, use maximum x anyways
        smallcell = true;
        tcell.fr = rectangle_fraction (tcell.ns, tcell.alphas, 
                    (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5});
        itr = rsolver_brent (&newx, xmin, xmax, &tcell, fit_volume_error, 
                             tolerance = 1e-14, maxitr = maxitr);
    }
    #else
    if (itr == -1)
        return oldx;
    #endif
    
    #if PRINTA
    fprintf(stderr, "AREA SOLVER done after %d itrs: oldx=%g newx=%g small cell = %d\n", 
        itr, oldx, newx, smallcell);
    #endif

    (void) smallcell;
    if (itr == maxitr || (itr == 0 && newx == HUGE)) {
        fprintf(stderr, "WARNING: area root solver didn't converge after %d iterations\n", itr);
        return oldx;
    }

    return newx;
}


double get_real_error (const void* data, double alpha) 
{
    tripoint tcell = (*(tripoint*)data);

    static const coord lhs = {-0.5,-0.5,-0.5}, rhs = {0.5,0.5,0.5};

    double fa = rectangle_fraction (tcell.nf, alpha, lhs, rhs);
    double frcalc = tcell.s*immersed_fraction (fa, tcell.nf, alpha, tcell.ns, 
                                              tcell.alphas, lhs, rhs, print = 0);
    #if PRINTA
    #if 1
    double err = frcalc - clamp(tcell.fr, 0, tcell.s);
    fprintf(stderr, "alpha = %0.15g fa = %0.15g frcalc = %0.15g tcell.fr = %0.15g err=%g tcell.s=%g\n", 
        alpha, fa, frcalc, tcell.fr, err, tcell.s);
    #endif
    #endif

    return frcalc - clamp(tcell.fr, 0, tcell.s);
}


extern scalar f;


double ghost_alpha (const tripoint tcell, double alphaMin, double alphaMax, int * numitr = NULL)
{
    static const coord lhs = {-0.5,-0.5,-0.5}, rhs = {0.5,0.5,0.5};
#if dimension == 2
    coord ip[2];
    int dp = facets (tcell.ns, tcell.alphas, ip);
#else
    coord ip[12];
    int dp = facets (tcell.ns, tcell.alphas, ip, 1);
#endif
    double alphaArray[4];
    double frArray[4];

    for (int i = 0; i < dp; ++i) {
        alphaArray[i] = 0;
        foreach_dimension()
            alphaArray[i] += tcell.nf.x*ip[i].x;
        double fi = plane_volume (tcell.nf, alphaArray[i]);
        frArray[i] = tcell.s*immersed_fraction (fi, tcell.nf, alphaArray[i], 
                                                tcell.ns, tcell.alphas, lhs, rhs);
    }

    for (int i = 0; i < dp; ++i) { // TODO: is there a better check than this?
        //if (approx_equal_double (frArray[i], tcell.fr, INT_TOL));
        if (frArray[i] >= tcell.fr - INT_TOL && frArray[i] <= tcell.fr + INT_TOL) {
        //if (frArray[i] >= tcell.fr - 1e-10 && frArray[i] <= tcell.fr + 1e-10) {
        #if PRINTA
            fprintf(stderr, "ghost_alpha: fr = %g alpha = %g fr0 = %g\n", 
                frArray[i], alphaArray[i], tcell.fr);
        #endif
            return alphaArray[i];
        }
    }

    // Fall back to basic root solver
    double alpha = 0;
    int itr = rsolver_brent (&alpha, alphaMin, alphaMax, &tcell, get_real_error, tolerance = 1e-12);

    if (numitr)
        *numitr = itr;

    return alpha;
}


trace
double immersed_alpha (double f, double ibm, coord nf, double alphaf, coord ns, double alphas,
                       double freal, double tolerance = BI_TOL, int * numitr = NULL)
{
    if ((freal <= VTOL || freal >= 1-VTOL) && !nf.x && !nf.y && !nf.z)
        return alphaf;

    //coord lhs = {-0.5,-0.5,-0.5}, rhs = {0.5,0.5,0.5}; // bottom-left & top-right points, resp.

    double alphaMin = -0.6, alphaMax = 0.6, alpha = 0; // are there better values?
 
    double f0 = clamp(f, 0., 1.);
    //double ibm0 = rectangle_fraction (ns, alphas, lhs, rhs);
    double ibm0 = ibm;
    freal = clamp(freal, 0., ibm0);
    const tripoint tcell = fill_tripoint (freal, nf, alphaf, ns, alphas, f0, ibm0);

    int maxitr = 50;

#if PRINTA
    fprintf(stderr, "|| A.S: nf={%g, %g, %g} ns={%g, %g, %g} alphas=%g\n",
        nf.x, nf.y, nf.z, ns.x, ns.y, ns.z, alphas);
#endif

    if (!nf.x && !nf.y && !nf.z)
        fprintf(stderr, "WARNING: all components of nf are 0 in alpha solver!\n");

#if 1
    // the cell is full or empty, but we want to change f to set C.A while conserving freal
      if ((nf.x || nf.y || nf.z) && (freal >= ibm0 - VTOL || freal <= 0 + VTOL)) {
        return ghost_alpha (tcell, alphaMin, alphaMax, numitr);
    }
#endif

    #if PRINTA
    fprintf (stderr, "|| A.S: f0=%0.15g alphaf=%0.15g freal=%0.15g ibm=%g\n", 
                      f0, alphaf, freal, ibm0);
    #endif

    //int itr = rsolver_bisection (&alpha, alphaMin, alphaMax, &tcell, get_real_error, maxitr = maxitr);
    int itr = rsolver_brent (&alpha, alphaMin, alphaMax, &tcell, get_real_error, 
                             maxitr = maxitr, tolerance = 1e-12);
    if (numitr)
        *numitr = itr;

    return alpha;
}


/*
this function fills fr with only the real volume fraction of f that lays 
outside of the solid immersed boundary, ibm.
*/

void real_fluid (scalar f, scalar fr, vector nf, scalar alphaf, vector ns, scalar alphas)
{
    foreach() {
        if (on_interface(ibm) && on_interface(f))
            fr[] = immersed_fraction (f[], (coord){nf.x[], nf.y[], nf.z[]}, alphaf[],
                                           (coord){ns.x[], ns.y[], ns.z[]}, alphas[],
                                           (coord){-0.5,-0.5,-0.5},
                                           (coord){ 0.5, 0.5, 0.5}) * ibm[];
        else
            fr[] = clamp(f[], 0., 1.) * ibm[];
    }
}


/*
interior cells are cells along the solid interface that
    1. are full before and after the advection step
    2. do not contribute to the contact angle imposition
    3. are exempt from the special three-phase advection treatment

TODO: make a formal definition of what an interior cell is (full cells along the solid interface?
      any full cell that would not contain an extrapolated interface for CAs?)
*/

bool is_interior_cell (Point point, scalar ibm, scalar c, scalar cr)
{
    #if PRINTA && 0
    fprintf(stderr, "checking if (%g, %g) is an interior cell ibm=%g cr=%g c=%g \n",
        x, y, ibm[], cr[], c[]);
    #endif

    //if ((ibm[] == 1 && cr[] == ibm[]) || cr[] == 0)
    //    return false;

    if (ibm[] <= 0 || ibm[] >= 1 || cr[] == 0)
        return false;

    //coord nf = interface_normal (point, c);
    //coord ns = interface_normal (point, ibm);

    //if (on_interface(ibm) && cr[] >= ibm[] - 1e-3 && !is_triple_point (point, nf, ns))
    //    return true;

    #if 0
    if (on_interface(ibm) && cr[] >= ibm[] - 1e-10)
        return true;

    #else
    foreach_neighbor() { // was (1)

        #if PRINTA && 0
        fprintf (stderr, "  || (%g, %g) ibm=%g cr=%g c=%g |cr - ibm| = %g\n",
            x, y, ibm[], cr[], c[], fabs(cr[] - ibm[]));
        #endif

    // TODO: if cr[] != ibm[], try seeing if ibm[] == 1 instead (full fluid cell)
        if (ibm[] && fabs(cr[] - ibm[]) > INT_TOL) { // CHANGED FROM 1e-3 TEMPORARILY
            return false;
        }
    }

    return true;
    #endif
    //return false;
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

double avgitr = 0;
int g_count = 0;

trace
void immersed_reconstruction (scalar c, scalar cr, vector nf, scalar alphaf, 
                              vector ns, scalar alphas, scalar id, int * fix = NULL)
{
    double error_sum = 0;
    int itrsum = 0, count = 0;
    foreach(reduction(+:error_sum) reduction(+:itrsum) reduction(+:count)) {
        c[] = clamp(c[], 0, 1);
        if (on_interface(ibm) && c[]) {

            #if PRINTA
            fprintf(stderr, "(%g, %g, %g) cr[] = %g c[] = %g ibm[] = %g\n", x, y, z, cr[], c[], ibm[]);
            #endif

            //if (cr[] < 1e-14 && c[] < 1e-14) {
            if (cr[] < 1e-12) {
                cr[] = c[] = 0;
                //id[] = -1;
                continue;
            }
            else if (approx_equal_double (cr[], ibm[]) && c[] >= 1.) {
                cr[] = ibm[];
                //id[] = -1;
                continue;
            }
            else if (is_interior_cell(point, ibm, c, cr)) {
                cr[] = ibm[];
                c[] = 1;
                continue;
            }

            cr[] = clamp(cr[], 0, ibm[]);

            coord nsolid = {ns.x[], ns.y[], ns.z[]}, nfluid = {nf.x[], nf.y[], nf.z[]};

            if (!nfluid.x && !nfluid.y && !nfluid.z) {
                c[] = 0.9999;
                nfluid = interface_normal (point, c);
            }

#if PRINTA
            fprintf(stderr, "nf={%g, %g, %g} af=%g ns={%g, %g, %g} as=%g c=%g cr=%g\n",
                nfluid.x, nfluid.y, nfluid.z, alphaf[], nsolid.x, nsolid.y, nsolid.z, 
                alphas[], c[], cr[]);
#endif

            double freal = immersed_fraction (c[], nfluid, alphaf[], nsolid, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0)*ibm[];
            //freal = clamp (freal, 0, ibm[]);

#if 0
            if (is_interior_cell (point, ibm, c, cr)) {
                c[] = 1;

                #if PRINTA
                fprintf(stderr, "(%g, %g) cell is interior cell c=%g cr=%g\n",
                    x, y, c[], cr[]);
                #endif

                continue;
            }
#endif
            error_sum += fabs(cr[] - freal) * pow(Delta, dimension);

            if (approx_equal_double (cr[], freal, 1e-10)) {
                //id[] = -1;
                continue;
            }

            int iter = 0;
            double alpha = immersed_alpha (c[], ibm[], nfluid, alphaf[], nsolid, alphas[], cr[], numitr = &iter);

            #if PRINTA
            double c0 = c[];
            #endif

            c[] = plane_volume (nfluid, alpha);
            alphaf[] = alpha;

            //id[] = 1;

            #if PRINTA
            fprintf(stderr, "(%g, %g, %g) c[]_before = %0.15f, c[]_after = %0.15f | cr[] = %0.15f crcalc =%0.15f\n",
                           x, y, z, c0, c[], cr[], freal);
            #endif

            itrsum += iter;
            count++;
        }
        if (cr[] < 1e-12 && ibm[] == 1) {
            cr[] = c[] = 0;
            //id[] = -1;
            continue;
        }
        //else
        //    id[] = -1;
    }

    boundary({id});

#if PRINTA
    fprintf(stderr, "error_sum = %g average # of iteration = %g count = %d\n", error_sum, ((double)itrsum)/(count+SEPS), count);
#endif

    avgitr = ((double)itrsum)/(count+SEPS);
    g_count = count;

    if (fix)
        *fix = error_sum > 1e-10? 1: 0;
}



trace
double immersed_alpha_temp (double f, double ibm, coord nf, double alphaf, coord ns, double alphas,
                            double freal, double tolerance = BI_TOL, int * numitr = NULL)
{
    if ((freal <= VTOL || freal >= 1-VTOL) && !nf.x && !nf.y && !nf.z)
        return alphaf;

    //coord lhs = {-0.5,-0.5,-0.5}, rhs = {0.5,0.5,0.5}; // bottom-left & top-right points, resp.

    double alphaMin = -0.6, alphaMax = 0.6, alpha = 0; // are there better values?
 
    double f0 = clamp(f, 0., 1.);
    double ibm0 = ibm;
    freal = clamp(freal, 0., ibm0);
    const tripoint tcell = fill_tripoint (freal, nf, alphaf, ns, alphas, f0, ibm0);

    int maxitr = 40;

#if PRINTA
    fprintf(stderr, "|| A.S: nf={%g, %g, %g} ns={%g, %g, %g} alphas=%g\n",
        nf.x, nf.y, nf.z, ns.x, ns.y, ns.z, alphas);
#endif

    if (!nf.x && !nf.y && !nf.z)
        fprintf(stderr, "WARNING: all components of nf are 0 in alpha solver!\n");

#if 1
    // the cell is full or empty, but we want to change f to set C.A while conserving freal
      //if ((nf.x || nf.y || nf.z) && (freal >= ibm0 - INT_TOL || freal <= 0 + INT_TOL)) {
      if ((nf.x || nf.y || nf.z) && (freal <= 0 + VTOL || freal >= ibm0 - VTOL)) {
        return ghost_alpha (tcell, alphaMin, alphaMax, numitr);
    }
#endif

    #if PRINTA
    fprintf (stderr, "|| A.S: f0=%0.15g alphaf=%0.15g freal=%0.15g ibm=%g\n", 
                      f0, alphaf, freal, ibm0);
    #endif

    //int itr = rsolver_bisection (&alpha, alphaMin, alphaMax, &tcell, get_real_error, maxitr = maxitr);
    int itr = rsolver_brent (&alpha, alphaMin, alphaMax, &tcell, get_real_error, maxitr = maxitr);
    if (numitr)
        *numitr = itr;

    return alpha;
}


trace
void immersed_reconstruction_temp (scalar c, scalar cr, vector nf, scalar alphaf, 
                                   vector ns, scalar alphas, int * fix = NULL)
{
    #if PRINTA
    fprintf(stderr, "TEMP RECONSTRUCTION!\n");
    #endif

    double error_sum = 0;
    int itrsum = 0, count = 0;
    foreach(reduction(+:error_sum) reduction(+:itrsum) reduction(+:count)) {
        c[] = clamp(c[], 0, 1);
        if (on_interface(ibm) && c[]) {

            #if PRINTA
            fprintf(stderr, "(%g, %g, %g) cr[] = %g c[] = %g ibm[] = %g\n", x, y, z, cr[], c[], ibm[]);
            #endif

#if 0
            if (cr[] < 1e-14 && c[] < 1e-14) {
                cr[] = c[] = 0;
                continue;
            }
#endif            
            if (approx_equal_double (cr[], ibm[]) && c[] >= 1.) {
                cr[] = ibm[];
                continue;
            }

            cr[] = clamp(cr[], 0, ibm[]);

            coord nsolid = {ns.x[], ns.y[], ns.z[]}, nfluid = {nf.x[], nf.y[], nf.z[]};

            if (!nfluid.x && !nfluid.y && !nfluid.z) {
                c[] = 0.9999;
                nfluid = interface_normal (point, c);
            }

#if PRINTA
            fprintf(stderr, "nf={%g, %g, %g} af=%g ns={%g, %g, %g} as=%g c=%g cr=%g\n",
                nfluid.x, nfluid.y, nfluid.z, alphaf[], nsolid.x, nsolid.y, nsolid.z, 
                alphas[], c[], cr[]);
#endif

            double freal = immersed_fraction (c[], nfluid, alphaf[], nsolid, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0)*ibm[];
            
            error_sum += fabs(cr[] - freal) * pow(Delta, dimension);

           if (approx_equal_double (cr[], freal, 1e-10)) {
                continue;
           }

            int iter = 0;
            double alpha = immersed_alpha_temp (c[], ibm[], nfluid, alphaf[], nsolid, alphas[], cr[], numitr = &iter);

            #if PRINTA
            double c0 = c[];
            #endif

            c[] = plane_volume (nfluid, alpha);
            //alphaf[] = alpha;


            #if PRINTA
            fprintf(stderr, "(%g, %g, %g) c[]_before = %0.15f, c[]_after = %0.15f | cr[] = %0.15f crcalc =%0.15f\n",
                           x, y, z, c0, c[], cr[], freal);
            #endif

            itrsum += iter;
            count++;
        }
    }

#if PRINTA
    fprintf(stderr, "error_sum = %g average # of iteration = %g count = %d\n", error_sum, ((double)itrsum)/(count+SEPS), count);
#endif

    avgitr = ((double)itrsum)/(count+SEPS);
    g_count = count;

    if (fix)
        *fix = error_sum > 1e-10? 1: 0;
}


/*
real_volume calculates the volume of the portion of f which lays outside of
the immersed boundary
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
            volume += immersed_fraction (f[], (coord){nf.x[], nf.y[], nf.z[]}, alphaf[],
                                              (coord){ns.x[], ns.y[], ns.z[]}, alphas[],
                                              (coord){-0.5, -0.5, -0.5},
                                              (coord){0.5, 0.5, 0.5}, 0) * dv2();
        else
            volume += f[]*dv2();
    }

    return volume;
}

#if CA

double get_contact_angle (scalar f, scalar ibm)
{
    vector nf[], ns[];
    scalar alphaf[], alphas[];

    reconstruction (f, nf, alphaf);
    reconstruction (ibm, ns, alphas);

    scalar fr_temp[];
    real_fluid (f, fr_temp, nf, alphaf, ns, alphas);

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
#endif 


/*
redistribute_volume calculates the volume loss caused by imprecise movement of the
solid boundary (ibm), which itself is caused by us using different methods
to advect the liquid and solid interface: VOF and level-set reinitialization -> VOF, resp.

The weights used for redistribution is based on the number of interfacial cells. Meaning
each cell will get the same amount of volume.

TODO: improve weighting function by maybe including distance from interface?
*/

double redistribute_volume (scalar c, scalar cr, const scalar ibm)
{
    // 1. Calculate the volume error
    //    and count the number of interfacial cells that contain no solid fragment
    double verror = 0;
    int icells = 0;
    foreach(reduction (+:verror) reduction (+:icells)) {
        if (cr[] > ibm[]) { // cell is too full
            verror += (cr[] - ibm[])*dv();
            cr[] = ibm[];
            c[] = 1;
        }
        else if (cr[] < 0) { // cell is too empty
            verror += cr[]*dv();
            c[] = cr[] = 0;
        }
        #if MOVING
        else if (on_interface(ibm) && cr[] < ibm[] && // moving interface error
            is_interior_cell (point, ibm, c, cr)) {
            verror += cr[] - ibm[]*dv();
            cr[] = ibm[];
        }
        #endif
        if (on_interface(c) && ibm[] >= 1)
            icells++;
    }
   
    if (verror == 0) return verror;

    boundary ({cr});

    // 2. Redistribute volume error to interfacial cells
    vector id[];

    int count = 0, overfill = 0;

    foreach(reduction(max:overfill) reduction(+:count)) {
        if (on_interface(c) && ibm[] >= 1) {

            // cell is too full to take the additional volume, so give it to neighbors
            if (cr[] + (verror/(icells*dv())) > ibm[]) {
                overfill = 1;

                count++;

                // check for left/right neighbors
                bool done = false;
                for (int i = -1; i <= 1; i += 2)
                    if (cr[i] + 2*(verror/(icells*dv())) < ibm[i] && ibm[i]) // 2* because cell may be a recipient for multiple cells (temp fix)
                        id.x[] = i, done = true;

                // then check for top/bottom neighbors
                if (!done)
                   for (int j = -1; j <= 1; j += 2)
                       if (cr[0,j] + 2*(verror/(pow(Delta,dimension)*cm[0,j]*icells)) < ibm[0,j] && ibm[0,j])
                           id.y[] = j, done = true;

#if dimension == 3
                if (!done)
                   for (int k = -1; k <= 1; k += 2)
                       if (cr[0,0,k] + 2*(verror/(icells*dv())) < ibm[0,0,k] && ibm[0,0,k])
                           id.z[] = k, done = true;
#endif
                
                if (!done)
                    fprintf (stderr, "WARNING: could not fill cell with volume error!\n");
            }
            else {

                cr[] += verror/(icells*dv());
                c[]  += verror/(icells*dv());
            }
        }
        else {
            foreach_dimension()
                id.x[] = 0;       
        }
    }

    boundary ({id});


    int fixed = overfill? 0: 1;
    if (overfill)
        foreach(reduction(max:fixed)) {
            for (int i = -1; i <= 1; i += 2)
                foreach_dimension()
                    if (id.x[i] == -i) {
                        //double cr0 = cr[];
                        cr[] += verror/(icells*pow(Delta,dimension)*cm[i]);
                        c[]  += verror/(icells*pow(Delta,dimension)*cm[i]);
                        fixed = 1;
                        #if 0
                        fprintf(stderr, "OVERFILL FIX:(%g, %g) cr0=%g cr[]=%g id=%g i=%d\n",
                            x, y, cr0, cr[], id.x[i], i);
                        #endif
                        cr[] = clamp (cr[], 0, ibm[]);
                        c[] = clamp (c[], 0, 1);
                    }
        }
#if PRINTA
    fprintf (stderr, "verror=%g icells=%d overfill=%d fixed=%d count=%d\n",
        verror, icells, overfill, fixed, count);
#endif
    boundary ({cr, c});
    (void) fixed; // to prevent unused variable warning

    return verror;
}

