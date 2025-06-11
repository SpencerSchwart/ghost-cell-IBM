/*
###### TWO PHASE FUNCTIONS FOR IBM ######
*/

#define PRINTA 0

#include "fractions.h"

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


/*
sort_clockwise sorts a list of coordinates, provided in cf w/nump points, in
clockwise order (or counter-clockwise if y-advection).
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
    qsort_r (cf, nump, sizeof(coord), is_behind, &pc);
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

    double minx = -0.5, maxx = 0, newx; 
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
    #if PRINTA
    fprintf(stderr, "AREA SOLVER done: oldx=%g newx=%g\n", oldx, newx);
    #endif
    if (itr == maxitr)
        fprintf(stderr, "WARNING: area root solver didn't converge after %d iteration wtih an error of %g\n",
                        itr, error);
    return newx;
}


/*
immersed_area calculates the area of the real fluid part of a cell given a 
bounding box (defined in lhs and rhs)
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
    #if 0 // temporarily off
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
    sort_clockwise (nump, cf, print);

    // 6. use the shoelace formula to find the area
    double area = polygon_area (nump, cf);

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


int bisection_solver (double* a, double amin, double amax, tripoint tcell, 
                      double tolerance = 1e-9, int maxitr = 40, int print = 0)
{
    coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}; // bottom-left & top-right points, resp.
   
    double ibm0 = tcell.s == 0? rectangle_fraction (tcell.ns, tcell.alphas, lhs, rhs): tcell.s;

    double cmax = rectangle_fraction (tcell.nf, amax, lhs, rhs);
    double crmax = immersed_fraction (cmax, tcell.nf, amax, tcell.ns, tcell.alphas, lhs, rhs)*ibm0;
    double errmax = crmax - tcell.fr;
    
    #if PRINTA
    fprintf(stderr, "|| A.S: cmax=%g crmax=%g errmax=%g\n", cmax, crmax, errmax);
    #endif

    if (fabs(errmax) < tolerance) { // problematic for receding, hydrophobic CAs
        *a = amax;
        return 0;
    }

    double alpha = 0, error = HUGE;
    int itr = 0;
    while (fabs(error) > tolerance && itr < maxitr) {
        alpha = (amin + amax)/2.;   // bisection
        double fa = rectangle_fraction (tcell.nf, alpha, lhs, rhs);
        double fcalc = immersed_fraction (fa, tcell.nf, alpha, tcell.ns, tcell.alphas, lhs, rhs, print)*ibm0;
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
        
    if (itr == maxitr) {
        fprintf(stderr, "WARNING: alpha  solver does not converge after"
                        " maximum iteration (%d), error = %g\n", maxitr, error);
    #if 0
        fprintf(stderr, "\n");
        for (double a = -0.75; a <= 0.75; a += 0.01) {
            double fa = rectangle_fraction (tcell.nf, a, lhs, rhs);
            double fcalc = immersed_fraction (fa, tcell.nf, a, tcell.ns, tcell.alphas, lhs, rhs)*ibm0;
            double error1 = fcalc - tcell.fr;
            fprintf (stderr, "%g %g %g %g\n", a, fa, fcalc, error1);
        }
        fprintf(stderr, "\n");

    #endif
    }
    *a = alpha;
    return itr;
}

extern scalar f;

double immersed_line_alpha (Point point, coord nf, double alphaf, coord ns, double alphas,
                            double freal, double tolerance = 1e-9)
{
    if ((freal <= INT_TOL || freal >= 1-INT_TOL) && !nf.x && !nf.y)
        return alphaf;

    coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}; // bottom-left & top-right points, resp.

    double alphaMin = -1, alphaMax = 1; // are there better values?
 
    double error = HUGE, alpha = 0;
    
    double f0 = clamp(f[], 0., 1.);
    double ibm0 = rectangle_fraction (ns, alphas, lhs, rhs);
    freal = clamp(freal, 0., ibm0);
    tripoint tcell = fill_tripoint (freal, nf, alphaf, ns, alphas, f0, ibm0);

    int maxitr = 40;

    // the cell is full or empty, but we want to change f to set C.A while conserving freal
      if ((nf.x || nf.y) && (freal >= ibm0 - INT_TOL || freal <= 0 + INT_TOL)) {

        coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}, ip[2]; // intersecting point
        int bp = boundary_points (ns, alphas, lhs, rhs, ip);
        double alpha0 = nf.x*ip[0].x + nf.y*ip[0].y;
        double alpha1 = nf.x*ip[1].x + nf.y*ip[1].y;

        double f0 = plane_volume (nf, alpha0);
        double f1 = plane_volume (nf, alpha1);

        double freal0 = ibm0*immersed_fraction (f0, nf, alpha0, ns, alphas, lhs, rhs);
        double freal1 = ibm0*immersed_fraction (f1, nf, alpha1, ns, alphas, lhs, rhs);
        if (freal0 >= freal - INT_TOL && freal0 <= freal + INT_TOL)
            return alpha0;
        else if (freal1 >= freal - INT_TOL && freal1 <= freal + INT_TOL)
            return alpha1;
        else {
            bisection_solver (&alpha, -0.75, 0.75, tcell);
            return alpha;
        }
    }
    alphaMin = -0.75, alphaMax = 0.75;

    #if PRINTA
    fprintf (stderr, "|| A.S: f0=%0.15g alphaf=%0.15g freal=%0.15g ibm=%g\n", 
                      f0, alphaf, freal, ibm0);
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
                                                 (coord){0.5, 0.5, 0.5}, 0) * ibm[];
        else
            fr[] = clamp(f[], 0., 1.) * ibm[];
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
        c[] = clamp(c[], 0, 1);
        if (on_interface(ibm) && c[]) {

            cr[] = clamp(cr[], 0, ibm[]);
            
            #if PRINTA
            fprintf(stderr, "cr[] = %g\n", cr[]);
            #endif

            coord nsolid = {ns.x[], ns.y[]}, nfluid = {nf.x[], nf.y[]};

            double freal = immersed_fraction (c[], nfluid, alphaf[], nsolid, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0)*ibm[];
            freal = clamp (freal, 0, ibm[]);
            
            if (cr[] < freal + INT_TOL && cr[] > freal - INT_TOL)
                continue;

            double alpha = immersed_line_alpha (point, nfluid, alphaf[], nsolid, alphas[], cr[]);
            double c0 = c[];
            c[] = plane_volume (nfluid, alpha);
            alphaf[] = alpha;

            #if PRINTA
            fprintf(stderr, "(%g, %g) c[]_before = %0.15f, c[]_after = %0.15f | cr[] = %0.15f crcalc =%0.15f\n",
                           x, y, c0, c[], cr[], freal);
            #endif
        }
    }
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
            volume += immersed_fraction (f[], (coord){nf.x[], nf.y[]}, alphaf[],
                                              (coord){ns.x[], ns.y[]}, alphas[],
                                              (coord){-0.5, -0.5, -0.5},
                                              (coord){0.5, 0.5, 0.5}, 0) * sq(Delta)*ibm[];
        else
            volume += f[]*sq(Delta)*ibm[];
    }

    return volume;
}

#if CA
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
#endif 

/*
real_interfacial checks to see if the given point should be used for extrapolation
while setting the CAs
*/

static inline bool real_interfacial (Point point, scalar cr, scalar c)
{
    if (cr[] >= ibm[] - INT_TOL && ibm[]) {
        for (int i = -1; i <= 1; i += 2)
            foreach_dimension() {
                if (cr[i] <= 0 && ibm[i])
                    return true;
                else if (cr[i] < ibm[i] - INT_TOL && cr[i] > INT_TOL && 
                         ibm[] > 0 && ibm[] < 1 && c[] > INT_TOL && c[] < 1-INT_TOL)
                    return true;
             }
    }
    else if (cr[] <= INT_TOL && ibm[]) {
        for (int i = -1; i <= 1; i += 2)
            foreach_dimension()
                if (cr[i] >= ibm[i] && ibm[i])
                    return true;
                else if (cr[i] < ibm[i] - INT_TOL && cr[i] > INT_TOL && 
                         ibm[] > 0 && ibm[] < 1 && c[] > INT_TOL && c[] < 1-INT_TOL)
                    return true;
    }
    else if (cr[] > INT_TOL && cr[] < ibm[] - INT_TOL)
        return true;
    return false;
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

