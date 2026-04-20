#define CA 1

#include "ibm-gcm-vof.h"

#include "mls.h"

#undef INT_TOL
#define INT_TOL 1e-10

scalar contact_angle[];

extern scalar f;
extern scalar ch;

extern double mu1; // liquid viscosity

attribute {
    struct cox_model {
        double ls, lm;           // microscopic and intermediate length scales (in meters)
    } cox_model;

    struct kistler_model {
        double ct;               // constant assumed to be about 72rad^3
    } kistler_model;

    struct shikh_model {
        double a2, a3, a4;       // phenomenological constants
    } shikh_model;

    struct wetting {
        double theta_s;          // static angle
        double theta_a, theta_r; // advancing and receding angles

        bool dynamic;            // to enable/disable dynamic CA model
        bool hysteresis_only;    // constant advancing/receding angles

        double (*model)(double, double, int); // returns the dynamic CA based on chosen model
    } wetting;
}

static inline double cox_ca_expression (double x)
{
    return cube(x)/9. - 0.00183985*pow(x, 4.5) + (1.845823E-6)*pow(x, 12.258487);
}

static inline double cox_ca_expression_inverse (double x)
{
    return pow(9*x, 1./3.) + 0.0727387*x - 0.0515388*sq(x) + 0.00341336*cube(x);
}

double cox_ca_model (double theta, double ca, int s)
{
    double val = cox_ca_expression(fabs(theta)) + s*ca*log(f.cox_model.lm/f.cox_model.ls);
    return cox_ca_expression_inverse(fabs(val));
}

static inline double kistler_ca_expression (double x)
{
    return acos(1 - 2*tanh(5.16*pow(x/(1 + 1.31*pow(x,0.99)), 0.706)));
}

static inline double kistler_ca_expression_inverse (double x)
{
    return cube(x)/f.kistler_model.ct;
}

static inline double receding_ca_model (double theta, double ca)
{
    double off = (cube(ca) - ca)/72.; // offset for a more smooth curve
    return pow(fabs(theta - 72*(ca - off)), 1./3.);
}

double kistler_ca_model (double theta, double ca, int s)
{
    double val;
    if (s > 0) // advancing
        val = kistler_ca_expression(ca + kistler_ca_expression_inverse(theta));
    else // receding
        val = receding_ca_model(theta, ca);

    return val;
}

double shikh_ca_model (double theta, double ca, int s)
{
    double val;
    if (s > 0) { // advancing 
        double a2 = f.shikh_model.a2;
        double a1 = 1 + (1 - a2)*(cos(theta) - f.shikh_model.a4);
        double u0 = (sin(theta) - theta*cos(theta))/(sin(theta)*cos(theta) - theta);
        double u1 = f.shikh_model.a3*ca;

        val = acos(clamp(cos(theta) - (2*u1*(a1 + a2*u0))/((1 - a2)*(sqrt(a1 + sq(u1)) + u1)), -1, 1));
    } else
        val = receding_ca_model(theta, ca);

    return val;
}

/**
Initialize the field storing the contact angles with the static contact angle first.*/
event defaults (i = 0)
{
    foreach()
        contact_angle[] = f.wetting.theta_s*pi/180.;

    contact_angle.refine = refine_injection;
    contact_angle.prolongation = refine_injection;
    boundary({contact_angle});

    if (!f.wetting.model) { // default if model not set
        f.wetting.model = cox_ca_model;
        f.cox_model.ls = 1e-9;
        f.cox_model.lm = 10e-6;
    }
}

/**
Basic 4-Point spreading funciton from Peskin.*/
double phi_func (double x, double h) 
{
    double r = x / h;
    double phi;

    if (fabs(r) <= 1)
      phi = (1./8.) * (3 - (2 * fabs(r)) + sqrt (1 + (4 * fabs(r)) - (4 * sq(r))));
    else if (fabs(r) > 1 && fabs(r) <= 2)
      phi = (1./8.) * (5 - (2 * fabs(r)) - sqrt (-7 + (12 * fabs(r)) - ( 4 * sq(r))));
    else
      phi = 0;
    return phi;
}

/**
delta_interpolation calculates the weight of a point (cpoint) with respect to a point (ipoint)
of which we want to interpolate to. This uses the dirac delta interpolation method.*/
double delta_interpolation (coord cpoint, coord ipoint, double dv, double delta)
{
    double dphi = 1.;
    foreach_dimension()
        dphi *= phi_func(cpoint.x - ipoint.x, delta);
    return dphi/dv;
}

/**
normal_contact returns the properly-oriented normal of an interface touching the 
immersed boundary. ns is the normal to the immersed solid boundary and nf is the
normal to the VOF interface not taking into account the contact angle boundary
condition. angle is the imposed static contact angle. */

static inline coord normal_contact (coord ns, coord nf, double angle)
{
    coord n;
#if dimension == 2
    if (- ns.x * nf.y + ns.y * nf.x >= 0) { // 2D cross product
        n.x = - ns.x * cos(angle) + ns.y * sin(angle);
        n.y = - ns.x * sin(angle) - ns.y * cos(angle);
    }
    else {
        n.x = - ns.x * cos(angle) - ns.y * sin(angle);
        n.y =   ns.x * sin(angle) - ns.y * cos(angle);
    }
#else // dimension == 3
#if 1
    foreach_dimension()  //Axis angle
      if (ns.x != 0) ns.x = -ns.x ; //We take the normal pointing toward the fluid

    double gamma = atan2(ns.z,ns.x); 
    double beta = asin(ns.y);  
    coord tau2 ={-sin(gamma), 0, cos(gamma)}, tau1;
    
    foreach_dimension() //cross product
        tau1.x = ns.y*tau2.z - ns.z*tau2.y; 

    double phi1 =0, phi2 =0;
    foreach_dimension(){
        phi1 += nf.x*tau1.x;
        phi2 += nf.x*tau2.x;
    }
    double phi = atan2(phi2,phi1); 
    
    //general configurations using Axis angle with spherical coordinates
    n.x = sin(angle)*cos(phi)*sin(beta)*cos(gamma) + cos(angle)*cos(beta)*cos(gamma) - 
          sin(angle)*sin(phi)*sin(gamma);
    n.y = cos(angle)*sin(beta) - sin(angle)*cos(phi)*cos(beta);
    n.z = sin(angle)*cos(phi)*sin(beta)*sin(gamma) + cos(angle)*cos(beta)*sin(gamma) + 
          sin(angle)*sin(phi)*cos(gamma);
#else
    foreach_dimension()
        if (ns.x != 0) ns.x = -ns.x;

    coord t;
    double cos_thetas = dot_product(nf, ns);
    foreach_dimension()
        t.x = nf.x - cos_thetas * ns.x;

    fprintf(stderr, "t={%g, %g, %g} ||t|| = %g", t.x, t.y, t.y, distance3D(t.x,t.y,t.z));

    normalize2(&t);

    foreach_dimension()
        n.x = cos(angle)*ns.x + sin(angle)*t.x;
#endif
#endif

    return n;
}


/**
Calculates the dynamic contact angle given ut (velocity tangent to the wall)
    - advancing: s =  1
    - receding:  s = -1 */

static inline double get_dynamic_ca (double theta, double sigma, double ut, int s)
{
    double capillary_num = mu1*fabs(ut)/sigma;
    double val = f.wetting.model(theta*pi/180., capillary_num, s);
    //fprintf(stderr, "val = %g\n", val);
    return val;
}

typedef struct hysteresis_error
{
    double c, s, alphas;
    coord pint, ns;
    int csign; // for orientating the interface correctly
} hysteresis_error;

double get_hysteresis_error (const void* data, double theta)
{
    const hysteresis_error tcell = (*(hysteresis_error*)data);
    coord ns = tcell.ns, nf;

    if (tcell.csign > 0) {
        nf.x = ns.x*cos(theta) - ns.y*sin(theta);
        nf.y = ns.x*sin(theta) + ns.y*cos(theta);
    }
    else {
        nf.x =  ns.x*cos(theta) + ns.y*sin(theta);
        nf.y = -ns.x*sin(theta) + ns.y*cos(theta);       
    }
    double alphaf = tcell.alphas + dot_product(nf, tcell.pint) - dot_product(ns, tcell.pint);

    static const coord lhs = {-0.5,-0.5,-0.5}, rhs = {0.5,0.5,0.5};

    double fa = rectangle_fraction (nf, alphaf, lhs, rhs);
    double frcalc = tcell.s*immersed_fraction (fa, nf, alphaf, tcell.ns, tcell.alphas, lhs, rhs);

    //fprintf(stderr, "%g nf={%g, %g} alphaf=%g f=%g fcalc=%g error=%g\n",
    //    180. - theta*180./pi, nf.x, nf.y, alphaf, fa, frcalc, frcalc - tcell.c);

    return frcalc - tcell.c;
}

double get_hysteresis_angle (double c, double s, coord pint, coord ns, double alphas, int csign)
{
    hysteresis_error tcell = {c, s, alphas, pint, ns, csign};

    double theta = 0;

    double theta_min = 0, theta_max = M_PI; // range

    rsolver_brent (&theta, theta_min, theta_max, &tcell, get_hysteresis_error, tolerance = 1e-10);
    return csign? M_PI - theta: theta;
}

scalar check[];

// ns is inward pointing normal of solid
int is_contact_cell (double c, double s, coord nf, coord ns, double alphas, double theta, coord * pint = NULL, double * alphaf = NULL)
{
    // calculate desired interface (normal, then intercept/alpha)
    normalize2(&ns);
    coord nc = normal_contact (ns, nf, theta);
    normalize_sum(&ns); normalize_sum(&nc);

    double alphac = immersed_alpha(c, s, nc, 0, ns, alphas, c);

    coord p;
    int points = interface_intersect(nc, alphac, ns, alphas, pint = &p);
    
    if (pint)
        *pint = p;

    if (alphaf)
        *alphaf = alphac;

    return points;
}

/**
Checks to see if a cell contains the contact point using the imposed contact angled boundary condition */
int is_contact_cell_bc (Point point, scalar c, scalar cs, face vector fs, coord * pint = NULL)
{
    // calculate desired interface (normal, then intercept/alpha)
    coord ns = facet_normal (point, cs, fs);
    double alphas = plane_alpha(cs[], ns);

    coord nf = interface_normal (point, c);

    normalize2(&ns);
    coord nc = normal_contact (ns, nf, contact_angle[]);
    normalize_sum(&ns); normalize_sum(&nc);

    double alphac = immersed_alpha(c[], cs[], nc, 0, ns, alphas, c[]);

    coord p;
    int points = interface_intersect(nc, alphac, ns, alphas, pint = &p);
    
    if (pint)
        *pint = p;

    return points;
}

/**
Calculates the midpoint of a three-phase cell by only considering the liquid line segment outside
the solid one (limited to 2D only). */
int contact_midpoint (coord nf, double alphaf, coord ns, double alphas, coord * p)
{
#if dimension == 2
    // face intersection points
    coord pc[2];
    int nump = facets(nf, alphaf, pc);

    if (nump < 2)
        return 0;
    
    // potential contact point with solid
    coord pint;
    if (interface_intersect (nf, alphaf, ns, alphas, pint = &pint)) {
        if (region_check2((plane){nf, alphaf}, (plane){ns, alphas}, pc[0]) > 0) // p[0] is inside solid
            pc[0] = pint;
        else
            pc[1] = pint;
    }

    foreach_dimension() 
        pint.x = 0.5*(pc[0].x + pc[1].x); 

    *p = pint;
#endif
    return 1;
}


void update_contact_angle (scalar c, scalar c0, scalar cs, vector ns, 
                           scalar alphas, vector u, scalar contact_angle)
{
    assert (dimension == 2);

    foreach() {
        check[] = 0;
        if (cs[] > 0 && cs[] < 1 && c[] > 0 && c[] < cs[] - INT_TOL) { // three-phase cell

#if 1
            bool receding = true;
            if (c0[] > 0 && c0[] < cs[] - INT_TOL) { // old three-phase cell

                // calculate the triple point coordinate of the original cell
                coord nf0 = interface_normal(point, c0), ns0 = {ns.x[], ns.y[]};

                normalize2(&ns0);
                coord nc0 = normal_contact (ns0, nf0, contact_angle[]);
                normalize_sum(&ns0); normalize_sum(&nc0);

                double alphac0 = immersed_alpha (c[], cs[], nc0, 0, ns0, alphas[], c0[]);

                coord nf = interface_normal(point, c);

                normalize2(&ns0);
                coord nc = normal_contact (ns0, nf, contact_angle[]);
                normalize_sum(&ns0); normalize_sum(&nc);

                double alphac = immersed_alpha (c[], cs[], nc0, 0, ns0, alphas[], c[]);

                coord pint0, pint1;
                bool was_tri = interface_intersect (nc0, alphac0, ns0, alphas[], pint = &pint0);
                bool is_new_tri = !was_tri && interface_intersect (nc, alphac, ns0, alphas[], pint = &pint1);

                if (is_new_tri) pint0 = pint1;

                //fprintf(stderr, "1 (%g, %g) %d %d %g %g %g n={%g,%g}\n", x, y, was_tri, is_new_tri, contact_angle[]*180/pi, c[], c0[], ns0.x, ns0.y);

                // check if the cell contains the contact point
                if (!was_tri && !is_new_tri) {
                    check[] = -1;
                } else { // contains the contact point
                    check[] = 5;
                    int csign = sign2(cross_product_2d(ns0, nc0)); // orient the interface correctly

                    double theta = get_hysteresis_angle(c[], cs[], pint0, ns0, alphas[], csign);

                    if (theta >= f.wetting.theta_r*pi/180.)
                        receding = false;

                    //fprintf(stderr, " 2 (%g, %g) c0=%g c=%g csign=%d theta = %g pint={%g, %g} CA=%g\n", 
                    //    x, y, c0[], c[],csign, theta*180./pi, pint0.x, pint0.y, contact_angle[]*180./pi);

                    /** 
                    If CA from hysteresis is within the adv. and rec. angles, then use it.
                    Otherwise, use the dynamic contact angle model to find CA. */
                    if (theta > f.wetting.theta_r*pi/180. && theta < f.wetting.theta_a*pi/180.) {
                        contact_angle[] = theta;
                        check[] = 1;
                        continue;
                    }
                    else if (f.wetting.hysteresis_only && theta < f.wetting.theta_r*pi/180.) {
                        contact_angle[] = f.wetting.theta_r*pi/180.;
                        check[] = 2;
                        continue;
                    }
                    else if (f.wetting.hysteresis_only && theta > f.wetting.theta_a*pi/180.) {
                        contact_angle[] = f.wetting.theta_a*pi/180.;
                        check[] = 3;
                        continue;
                    }
                }
            }
#endif
            if (!f.wetting.hysteresis_only && contact_angle[] != 0) {
                coord nf = interface_normal(point, c), ns0 = {ns.x[], ns.y[]}, ts = {-ns.y[], ns.x[]};

                normalize2(&ns0); normalize2(&ts);
                coord nc = normal_contact (ns0, nf, contact_angle[]);
                normalize_sum(&ns0); normalize_sum(&nc);

                double alphac = immersed_alpha (c[], cs[], nc, 0, ns0, alphas[], c[]);

                coord pint;
                if (!contact_midpoint (nc, alphac, ns0, alphas[], p = &pint))
                    continue;
                
                pint = lines_intersect(nc, alphac, ns0, alphas[]);
                coord pintg = local_to_global_coord(point, pint);
                //coord pintg0 = pintg;

                // shift interpolation point
                int s = receding? -1: 1, o = 1;
 #if 1              
                normalize2(&nc); normalize(&ns0);
                if (cross_product_2d(ns0, nc) > 0) // orient correctly
                    o = -1;
                
                pintg = (coord){pintg.x + o*ts.x*(1.5*Delta), pintg.y + o*ts.y*(1.5*Delta)}; // shift 1.5 cells into liquid
                pintg = (coord){pintg.x + -ns0.x*(0.5*Delta), pintg.y + -ns0.y*(0.5*Delta)}; // shift 0.5 cells above wall
                
                
                //fprintf(stderr, "3 pintg = {%g, %g} pintg0={%g, %g nf={%g, %g} nc={%g, %g} ac=%g ns={%g, %g} as=%g ts={%g, %g} ca=%g\n", 
                //    pintg.x, pintg.y, pintg0.x, pintg0.y, nf.x, nf.y, nc.x, nc.y, alphac, ns0.x, ns0.y, alphas[], ts.x, ts.y, contact_angle[]*180/pi);

                int count = 0;
                coord pstencil[25];
                double ustencil[25];
                foreach_neighbor() {
#if IBM
                    if ((ch[] >= 0.5 && cs[] == 1) || (u.x[] && ch[] && cs[] <= 0.5)) {
#elif EMBED
                    if ((ch[] >= 0.5 && u.x[])) {
#endif
                        pstencil[count] = (coord){x,y,z};
                        ustencil[count] = dot_product((coord){u.x[], u.y[]}, ts);
                        count++;
                    }
                }
                coord * pcells = malloc(count * sizeof(coord));
                double * ucells = malloc(count * sizeof(double));

                for (int i = 0; i < count; i++) {
                    pcells[i] = pstencil[i];
                    ucells[i] = ustencil[i];
                }

                if (count < 3)
                    continue;

                double ut = interpolate_mls(count, pcells, ucells, pintg, 3*Delta, Delta, weighting_func = gaussian_function);

                if (pcells) free(pcells);
                if (ucells) free (ucells);

                if (ut == 0)
                    continue;
#endif

                double theta = receding? f.wetting.theta_r: f.wetting.theta_a;
                contact_angle[] = get_dynamic_ca(theta, c.sigma, ut, s);
                
                //fprintf(stderr,"4 (%g, %g) %d %d %g %g %g\n", x, y, receding, count, ut, theta, contact_angle[]*180/pi);

                check[] = check[] == -1 || check[] == 0? -1: 4;
            }
        }
    }

    /**
    If a three-phase cell does not contain the triple point, we want to make sure it
    takes on the contact angle of the cell that does, assuming its nearby.

    If no cell contains the contact point nearby, then we assume its not pinned and the
    dynamic ca model is used. */
    #if 1
    foreach() {
        if (check[] == -1 || check[] == 4) {
            double theta = contact_angle[];
            //double theta0 = theta;
            bool pinned = false;
            foreach_neighbor() {
                if (check[] == 1 || check[] == 2 || check[] == 3) { // use CA of pinned cell
                    theta = contact_angle[];
                    pinned = true;
                }
            }

            if (check[] == -1 && !pinned && !f.wetting.hysteresis_only) {
                foreach_neighbor()
                    if (check[] == 4)
                        theta = contact_angle[];
            }

            contact_angle[] = theta;
            //fprintf(stderr, "fix (%g, %g) %g %d %g %g\n", x, y, check[], pinned, theta0, contact_angle[]);
        }
        if (check[] == 1) {
            double theta = contact_angle[];
            foreach_neighbor() { // check to see if CL is actually pinned or not
                if (check[] == 3 || check[] == 2) {
                    theta = contact_angle[];
                }
            }
            contact_angle[] = theta;
        }
    }
    #endif
    boundary({contact_angle});
}


/**
clean_fluid is used to remove any fluid inside the solid boundary that is
not necessary in enforcing the contact angle. */

void clean_fluid (scalar f, scalar cs)
{
    foreach() {
        if (cs[] == 0. && f[]) {
            int fluidNeighbors = 0;
            foreach_neighbor() {
                if (f[] > VTOL && cs[]) {
                    ++fluidNeighbors;
                    break;
                }
            }
            if (!fluidNeighbors) 
                f[] = 0;
        }
    }
}

void clean_fluid_real (scalar f, scalar fr, scalar cs)
{
    foreach() {
        if (cs[] == 0. && f[]) {
            int fluidNeighbors = 0;
            foreach_neighbor() {
                if (fr[] > VTOL && cs[]) {
                    ++fluidNeighbors;
                    break;
                }
            }
            if (!fluidNeighbors) 
                f[] = 0;
        }
    }
}

/** 
 */
bool contact_line_border (Point point, scalar c, scalar cs, double angle, double tolerance = 1e-2)
{

    int vol = angle >= pi/2.? 1: 0;

    for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++) {
            if ((i == 0 && j == 0) || abs(i) == abs(j))
                continue;
            else {
                if (cs[i,j] > 0 && cs[i,j] < 1 && approx_equal_double(c[i,j], vol, tolerance))
                    return true;
            }
        }

#if 0
    /** check left-right */
    for(int i = -1; i <= 1; i += 2)
       if (cs[i] > 0 && cs[i] < 1 && approx_equal_double(c[i], vol*cs[i], tolerance)) 
           return true;

    /** check top-bottom */
    for(int j = -1; j <= 1; j += 2)
       if (cs[0,j] > 0 && cs[0,j] < 1 && approx_equal_double(c[0,j], vol*cs[0,j], tolerance))
           return true;

#if dimension == 3
    /** check front-back */
    for(int k = -1; k <= 1; k += 2)
       if (cs[0,0,k] > 0 && cs[0,0,k] < 1 && approx_equal_double(c[0,0,k], vol*cs[0,0,k], tolerance))
           return true;
#endif
#endif
    return false;
}

static inline bool has_n_border (Point point, scalar c, double tolerance = INT_TOL)
{
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if ((i == 0 && j == 0) || (abs(i) == abs(j)))
                continue;
            else {
                if (cs[i,j] && c[i,j] <= tolerance)
                    return true;
            }
        }
    }
    return false;
}

static inline bool has_border (Point point, scalar c, double tolerance = INT_TOL)
{

    for (int i = -1; i <= 1; i++)
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0)
                continue;
            else {
                if (c[i,j] >= tolerance)
                    return true;
            }
        }


#if 0
    /** check left-right */
    for(int i = -1; i <= 1; i += 2)
       if (c[i] >= tolerance)
           return true;

    /** check top-bottom */
    for(int j = -1; j <= 1; j += 2)
       if (c[0,j] >= tolerance)
           return true;

#if dimension == 3
    /** check front-back */
    for(int k = -1; k <= 1; k += 2)
       if (c[0,0,k] >= tolerance)
           return true;
#endif
#endif

    return false;
}


#include "heights.h"

//#include "curvature.h"

/**
This function extends the normal reconstruction() function in fractions.h by
taking into account the contact angle on immersed boundaries. Given a volume fraction
field, f, it fills fields n and alpha with the each cells normal and intercept, resp. 

TODO: clean up function! e.g., better variable names and remove unnecessary variables*/

vector hg[];

vector hfg0[], hfg1[], hfg2[];
//scalar kappag0[];

void reconstruction_contact (scalar c, scalar cr, vector n, scalar alpha, 
                                  vector ns, scalar alphas, scalar inter, 
                                  scalar ghostInter, scalar extra)
{
    //heights(c, hfg0);

#if 0
    scalar c0[];
    foreach() {
        if (cs[] > 0 && cs[] < 1 && (cr[] <= 1e-6 || cr[] >= cs[] - 1e-6) &&
            !is_interior_cell(point, cs, cr) && level == depth())
        {
            int near = 0;
            foreach_neighbor() //TODO: this is probably not necessary
                if (c[] > 0 && cs[]) { near = 1; break; }
            ghostInter[] = near;
            //ghostInter[] = 0;
        }
        else
            ghostInter[] = 0;
    }
#endif
    foreach() {
        //inter[] = c[] > 0 && c[] < cs[];
        inter[] = cr[] > VTOL && cr[] < cs[] - INT_TOL;
#if 1
        if (inter[] && contact_angle[] && full_interfacial(point, cs)) {
           //&& (has_n_border(point, c, 0.))) {
            if (cs[] > 0 && cs[] < 1) {
                extra[] = (int)contact_line_border(point, c, cs, contact_angle[], 1e-2);
                //extra[] = 1;
            }
            else { // full cell
                bool caCell = false;
                foreach_neighbor() {
                    if (cs[] && cr[]/cs[] <= INT_TOL) {
                        caCell = true;
                        break;
                    }
                }
                extra[] = caCell;
            }
        }
        else
            extra[] = 0;
#else
        extra[] = inter[] && cs[] > 0 && cs[] < 1 && contact_angle[];

#if 1 // Make sure that extrapolation cells isn't an interior cell
        if (extra[]) {
            bool caCell = false;
            foreach_neighbor() {
                if (cs[] && cr[]/cs[] < 0.9) {
                  caCell = true;
                  break;
                }
            }
            extra[] = caCell;
        }
#endif
#endif
    }

    foreach() {
        if (extra[]) {
            coord ns0;
            double alphas0 = 0;
            if (cs[] == 1) {
                ns0 = full_interfacial_normal(point, cs, &alphas0);
            }
            else {
                ns0 = (coord){ns.x[], ns.y[], ns.z[]};
                alphas0 = alphas[];
            }

            coord nf = {n.x[], n.y[], n.z[]};
            if (!nf.x && !nf.y && !nf.z) {
                extra[] = 0;
                continue;
            }

            normalize (&nf);
            normalize2 (&ns0);
            coord nc = normal_contact (ns0, nf, contact_angle[]);
            normalize_sum(&nc);
            normalize_sum(&ns0);

            foreach_dimension()
                n.x[] = nc.x;

            alpha[] = immersed_alpha(c[], cs[], nc, alpha[], ns0, alphas0, cr[]);

            coord pint = lines_intersect (nc, alpha[], ns0, alphas0);

            //if (fabs(pint.x) >= 0.6 || fabs(pint.y) >= 0.6)
            //   extra[] = 0;

            c[] = max(plane_volume(nc, alpha[]), cr[]);

            //fprintf(stderr, "(%g, %g) c=%g nc={%g, %g} ac=%g ns={%g, %g} as=%g cr=%g\n",
            //    x, y, c[], n.x[], n.y[], alpha[], ns.x[], ns.y[], alphas[], cr[]);

            #if 0
            if (c[] <= 0) { // so ch and cr are consistent! prevents sudden removal of gginter
                cr[] = 0;
            }
            #endif
#if 0
                if (!is_interior_cell(point, cs, cr) && level == depth()) {
                    int near = 0;
                    foreach_neighbor() 
                        if (c[] > 0 && cs[]) { near = 1; break; }
                    ghostInter[] = near;
                }
            }
            inter[] = c[] > 0 && c[] < 1;
            extra[] = inter[] && cs[] > 0 && cs[] < 1 && !ghostInter[] && cr[] < cs[]-1e-6 && contact_angle[];
#endif            
        }
    }

    foreach() {
        ghostInter[] = 0;
        if (on_interface(cs) && !extra[] && level == depth()) {
        //if (on_interface(cs) && (cr[] <= VTOL || cr[] >= cs[] - VTOL) && !extra[] && 
        //    !is_interior_cell(point, cs, cr)) {
            bool isghost = false;
            foreach_neighbor() {
                if (extra[]) {
                    isghost = true;
                    break;
                }
            }
            if (isghost)
                ghostInter[] = 1;
        }
    }

    extra.dirty = true;
    boundary ({c, cr, n, alpha, inter, ghostInter, extra});
}

/**
This is the event that modifies the ghost cell's volume fraction, ghostf, at each timestep
necessary to impose the desired contact angle while conserving the real fluid, fr. */

scalar gginter[];
scalar gf00[], gf0[], gf1[];
scalar inter[];
scalar extra[];

trace
void set_contact_angle (scalar c, scalar cr0, const scalar cs,
                        vector nf, scalar alphaf, vector ns, scalar alphas) // clean up names, more consistent (c -> ch)
{
    scalar ghostInter[];

    foreach() { gf00[] = c[]; }

    reconstruction_contact(c, cr0, nf, alphaf, ns, alphas, inter, ghostInter, extra);
    scalar c0[], c1[];
    scalar_clone (c0, c);
    scalar_clone (c1, c);

    #if 1
    foreach() {
        gginter[] = ghostInter[];
        c0[] = c[];
        gf0[] = c0[];
    }
    boundary({c0});
    #endif

    //heights(c, hfg1);

    // 3. look for ghost and pure solid cells near triple point to enforce B.C
    foreach() {
        c0[] = c[];
        if (cs[] <= 0. && level == depth()) {  // pure solid cell
            double ghostf = 0., totalWeight = 0.;
            coord ghostCell = {x,y,z};

            // 3.a. Calculate weights
            foreach_neighbor() {
                if (extra[] && is_active(cell)) {
                    // why cr0? to mitigate problems for a very slow moving contact line and
                    // a new cell only gets a little bit of f (due to the low velocity/flux)
                    double cellWeight = cs[] == 1? cr0[] * (1 - cr0[]): cs[] * (1. - cs[]) * cr0[] * (cs[] - cr0[]);

                    double dc = sqrt(sq(ghostCell.x - x) + sq(ghostCell.y - y) + sq(ghostCell.z - z))/Delta;
                    cellWeight /= pow(dc,15);

                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf1 = {nf.x[], nf.y[], nf.z[]};

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                }
            }

            // 3.b. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.) 
                c[] = ghostf / (totalWeight + SEPS);
#if 0
            else if (contact_angle[] <= 0.55*pi) {
                bool check = false;
                foreach_neighbor(1) {
                    if (cs[] && cr0[] && !extra[])
                        check = true;
                    else if (cs[] && extra[] ) {
                        check = false;
                        break;
                    }
                }
                if (check)
                    c[] = 1;
                else {
                    foreach_neighbor() {
                        if (cs[] && cr0[] && !extra[])
                            check = true;
                        else if (cs[] && extra[]) {
                            check = false;
                            break;
                        }
                    }
                    if (check)
                        c[] = 1;
                }
            }
#endif            
        }
        if (ghostInter[])
        {
            double ghostf = 0., totalWeight = 0.;
            coord ghostCell = {x,y,z};

            double gr0 = gf0[];
            double fr0 = cr0[];
            double cs0 = cs[];

            // 3.c. extrapolate f from direct neighbors to get initial guess
            //int cond0 = cr0[] <= INT_TOL;
            coord ns0 = {ns.x[], ns.y[], ns.z[]}, nf1;
            double alphas0 = alphas[];

            foreach_neighbor() {
                if (extra[] && is_active(cell)) {
                    nf1 = (coord){nf.x[], nf.y[], nf.z[]};

#if 0
                    if ((contact_angle[] > 3.*pi/4. || contact_angle[] < pi/4.) &&
                        cond0 && dot_product_norm(nf1, ns0) >= 0)
                        continue;
#endif
                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }
                    double fa = rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                    double alpha = plane_alpha(fa, nf1);

                    normalize_sum(&ns0);
                    double fr = immersed_fraction(fa, nf1, alpha, ns0, alphas0, (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5});
#if 1
                    
                    //if (!approx_equal_double(gr0, fr) && !(contact_angle[] < pi/2. && fr0 <= VTOL)) {
                    //if (!approx_equal_double(gr0, fr) && fr < gr0) { //TODO: THIS IS BAD FOR GENERAL CASES IMPROVE
                    if (!approx_equal_double(gr0, fr)) { //TODO: THIS IS BAD FOR GENERAL CASES IMPROVE

                        const tripoint tcell = fill_tripoint (fr0, nf1, alpha, ns0, alphas0, fa, cs0);
                        double alphatmp = ghost_alpha(tcell, -0.6, 0.6);
                        fa = plane_volume(nf1, alphatmp);
                        //fprintf(stderr, "(%g, %g) (%g, %g) gr0=%g fr=%g nf={%g, %g} af=%g ns={%g, %g} as=%g cr=%g fa=%g\n",
                        //    ghostCell.x, ghostCell.y, x, y, gr0, fr, nf.x[], nf.y[], alpha, ns.x[], ns.y[], alphas[], cr0[], fa);
                    }
#endif
                    double cellWeight = cs[] == 1? cr0[] * (1 - cr0[]): cs[] * (1. - cs[]) * cr0[] * (cs[] - cr0[]);

                    double dc = sqrt(sq(ghostCell.x - x) + sq(ghostCell.y - y) + sq(ghostCell.z - z))/Delta;
                    cellWeight /= pow(dc,15);

                    totalWeight += cellWeight;

                    ghostf += cellWeight * fa;
                }
            }

            if (totalWeight > 0.) {
                c[] = max(ghostf / totalWeight, cr0[]);
                if (c[] > 1 - INT_TOL) c[] = 1; // is this good?
            }

            #if 0
            // 3.e otherwise, cell should not be used to enforce C.A.
            else if (ghostf <= 0 && cr0[] <= INT_TOL) {
                c[] = 0;
            }
            #endif
        }
        c1[] = c[];
    }

    //heights(c, hfg2);

    //curvature(c, kappag0, f.sigma, add = false);

#if 0
    scalar ctmp[];
    foreach() 
        ctmp[] = c[];
    boundary({ctmp,cr0,c});
#endif

#if 0
    /**
    TODO: is this necessary for ALL contact angles? or just super hydrophobic?*/
    foreach() {
        gf1[] = c[];
        // make full cells with fractional ch conserve cr
        if (ghostInter[] && cr0[] >= cs[] - 1e-6 && c[] != c0[] && c[] != 1.) {
            coord mf = interface_normal(point, c1);
            if (!mf.x && !mf.y && !mf.z) {
                double ctemp = c1[];
                c1[] = 0.5;
                mf = interface_normal(point, c1);
                c1[] = ctemp;
            }
            coord ns1 = {ns.x[], ns.y[], ns.z[]};
            double alpha = plane_alpha(c1[], mf);           
            double fr = immersed_fraction(c1[], mf, alpha, ns1, alphas[], 
                                         (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5})*cs[];
            
            if (!approx_equal_double(fr, cr0[])) {
                const tripoint tcell = fill_tripoint (cr0[], mf, alpha, ns1, alphas[], c[], cs[]);
                double alphatmp = ghost_alpha(tcell, -0.6, 0.6);
                c[] = plane_volume(mf, alphatmp);
            }
        }
        #if 1
        if (extra[]) {
            coord nf1 = interface_normal(point, c1);
            if (!nf1.x && !nf1.y && !nf1.z) {
              c[] = 0;
              continue;
            }
            coord ns1 = {ns.x[], ns.y[], ns.z[]};
            double alpha = immersed_alpha(c[], cs[], nf1, alphaf[], ns1, alphas[], cr0[]);
            c[] = plane_volume(nf1, alpha);
        }
        #endif
    }
#endif
    //clean_fluid (c, cs);
    boundary ({c, cr0});
}

