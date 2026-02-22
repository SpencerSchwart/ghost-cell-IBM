#define CA 1

#include "ibm-gcm-vof.h"

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
    return pow(theta - 72*(ca - off), 1./3.);
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

        val = acos(cos(theta) - (2*u1*(a1 + a2*u0))/((1 - a2)*(sqrt(a1 + sq(u1)) + u1)));
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
    if (- ns.x * nf.y + ns.y * nf.x > 0) { // 2D cross product
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

    //fprintf(stderr, "starting root solver\n");
    int itr = rsolver_brent (&theta, theta_min, theta_max, &tcell, get_hysteresis_error, tolerance = 1e-10);
    //fprintf(stderr, "rsolver done in %d iterations\n", itr);

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


void update_contact_angle (scalar c, scalar c0, scalar ibm, vector ns, 
                           scalar alphas, vector u, scalar contact_angle)
{
    assert (dimension == 2);

    foreach() {
        check[] = 0;
        if (ibm[] > 0 && ibm[] < 1 && c[] > 0 && c[] < ibm[] - INT_TOL) { // three-phase cell

#if 1
            bool receding = true;
            if (c0[] > 0 && c0[] < ibm[] - INT_TOL) { // old three-phase cell

                // calculate the triple point coordinate of the original cell
                coord nf0 = interface_normal(point, c0), ns0 = {ns.x[], ns.y[]};

                normalize2(&ns0);
                coord nc0 = normal_contact (ns0, nf0, contact_angle[]);
                normalize_sum(&ns0); normalize_sum(&nc0);

                double alphac0 = immersed_alpha (c[], ibm[], nc0, 0, ns0, alphas[], c0[]);

                coord nf = interface_normal(point, c);

                normalize2(&ns0);
                coord nc = normal_contact (ns0, nf, contact_angle[]);
                normalize_sum(&ns0); normalize_sum(&nc);

                double alphac = immersed_alpha (c[], ibm[], nc0, 0, ns0, alphas[], c[]);

                coord pint0, pint1;
                bool was_tri = interface_intersect (nc0, alphac0, ns0, alphas[], pint = &pint0);
                bool is_new_tri = !was_tri && interface_intersect (nc, alphac, ns0, alphas[], pint = &pint1);

                if (is_new_tri) pint0 = pint1;

                fprintf(stderr, "1 (%g, %g) %d %d %g %g %g n={%g,%g}\n", x, y, was_tri, is_new_tri, contact_angle[]*180/pi, c[], c0[], ns0.x, ns0.y);

                // check if the cell contains the contact point
                if (!was_tri && !is_new_tri) {
                    check[] = -1;
                } else { // contains the contact point
                    check[] = 5;
                    int csign = sign2(cross_product_2d(ns0, nc0)); // orient the interface correctly

                    double theta = get_hysteresis_angle(c[], ibm[], pint0, ns0, alphas[], csign);

                    if (theta >= f.wetting.theta_r*pi/180.)
                        receding = false;

                    fprintf(stderr, " 2 (%g, %g) c0=%g c=%g csign=%d theta = %g pint={%g, %g} CA=%g\n", 
                        x, y, c0[], c[],csign, theta*180./pi, pint0.x, pint0.y, contact_angle[]*180./pi);

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

                double alphac = immersed_alpha (c[], ibm[], nc, 0, ns0, alphas[], c[]);

                coord pint;
                if (!contact_midpoint (nc, alphac, ns0, alphas[], p = &pint))
                    continue;
                
                pint = lines_intersect(nc, alphac, ns0, alphas[]);
                coord pintg = local_to_global_coord(point, pint);
                coord pintg0 = pintg;

                // shift interpolation point
                int s = receding? -1: 1, o = 1;
 #if 1              
                normalize2(&nc); normalize(&ns0);
                if (cross_product_2d(ns0, nc) > 0) // orient correctly
                    o = -1;
                
                pintg = (coord){pintg.x + o*ts.x*(1.5*Delta), pintg.y + o*ts.y*(1.5*Delta)}; // shift 1.5 cells into liquid
                pintg = (coord){pintg.x + -ns0.x*(0.5*Delta), pintg.y + -ns0.y*(0.5*Delta)}; // shift 0.5 cells above wall
                
                
                fprintf(stderr, "3 pintg = {%g, %g} pintg0={%g, %g nf={%g, %g} nc={%g, %g} ac=%g ns={%g, %g} as=%g ts={%g, %g} ca=%g\n", 
                    pintg.x, pintg.y, pintg0.x, pintg0.y, nf.x, nf.y, nc.x, nc.y, alphac, ns0.x, ns0.y, alphas[], ts.x, ts.y, contact_angle[]*180/pi);
                   #endif
#if 0
                double weight_correction = 0;
                foreach_neighbor() {
                    if (ch[] >= 0.5 && (u.x[] || u.y[]))
                        weight_correction += delta_interpolation((coord){x,y,z}, pintg, dv(), Delta)*dv();
                }
                //fprintf(stderr, "weight_correction = %g\n", weight_correction);


                double ut_sum = 0, weight_sum = 0;
                coord ui = {0,0,0};

                foreach_neighbor() { // average velocity in interfacial cells to get contact line velocity

                    if (ch[] >= 0.5 && (u.x[] || u.y[])) {
                        double weight = delta_interpolation((coord){x,y,z}, pintg, dv(), Delta);
                        weight_sum += weight;

                        foreach_dimension()
                            ui.x += u.x[] * weight * dv() / weight_correction;

                        //fprintf(stderr, "%d (%g, %g) ut_sum=%g weight_sum=%g weight=%g ut=%g pint={%g,%g} pintg={%g,%g}\n", 
                        //    receding, x, y, ut_sum, weight_sum, weight * dv(), dot_product((coord){u.x[],u.y[],0}, ts), pint.x, pint.y, pintg.x, pintg.y);
                    }
                }
                double ut = -dot_product(ui, ts);
                //fprintf(stderr, "ut = %g\n", ut);
#else
                int count = 0;
                coord pstencil[25];
                double ustencil[25];
                foreach_neighbor() {
                    //double r = sqrt(sq(x - pintg.x) + sq(y - pintg.y));
                    if ((ch[] >= 0.5 && ibm[] == 1) || (u.x[] && ch[] && ibm[] <= 0.5)) {
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

                if (count == 0)
                    continue;

                double ut = interpolate_mls(count, pcells, ucells, pintg, 3*Delta, Delta, weighting_func = gaussian_function);

                if (pcells) free(pcells);
                if (ucells) free (ucells);
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
            double theta0 = theta;
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

#if 0
            if (!pinned && !f.wetting.hysteresis_only) {
                coord nf = interface_normal(point, c), ns0 = {ns.x[], ns.y[]}, ts = {-ns.y[], ns.x[]};

                normalize2(&ns0); normalize2(&ts);
                coord nc = normal_contact (ns0, nf, contact_angle[]);
                normalize_sum(&ns0); normalize_sum(&nc);

                int csign = sign2(cross_product_2d(ns0, nc)); // orient the interface correctly

                double ut_sum = 0;
                int count = 0;
                foreach_neighbor() { // average velocity in interfacial cells to get contact line velocity
                    if (ibm[] == 1 && c[] > INT_TOL && c[] < 1 - INT_TOL) {
                        ut_sum += dot_product ((coord){u.x[],u.y[],0}, ts);
                        count++;
                    }
                }

                if (count == 0) {
                    contact_angle[] = 0;
                    continue;
                }

                double ut = -ut_sum /((double)count);
                
                double theta = f.wetting.theta_a;
                int s;
                if (csign == sign(ut)) { 
                    theta = f.wetting.theta_r; // receding CL
                    s = -1;
                }
                else {
                    theta = f.wetting.theta_a; // advancing CL
                    s = 1;
                }
                contact_angle[] = get_dynamic_ca(theta, c.sigma, ut, s);
                //check[] = 4;
            }
#endif
        }
    }
    #endif
    boundary({contact_angle});
}


/**
clean_fluid is used to remove any fluid inside the solid boundary that is
not necessary in enforcing the contact angle. */

void clean_fluid (scalar f, scalar ibm)
{
    foreach() {
        if (ibm[] == 0. && f[]) {
            int fluidNeighbors = 0;
            foreach_neighbor() {
                if (f[] > VTOL && ibm[]) {
                    ++fluidNeighbors;
                    break;
                }
            }
            if (!fluidNeighbors) 
                f[] = 0;
        }
    }
}

void clean_fluid_real (scalar f, scalar fr, scalar ibm)
{
    foreach() {
        if (ibm[] == 0. && f[]) {
            int fluidNeighbors = 0;
            foreach_neighbor() {
                if (fr[] > VTOL && ibm[]) {
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
This function extends the normal reconstruction() function in fractions.h by
taking into account the contact angle on immersed boundaries. Given a volume fraction
field, f, it fills fields n and alpha with the each cells normal and intercept, resp. 

TODO: clean up function! e.g., better variable names and remove unnecessary variables*/

void reconstruction_contact (scalar c, scalar cr, vector n, scalar alpha, 
                                  vector ns, scalar alphas, scalar inter, 
                                  scalar ghostInter, scalar extra)
{
    scalar c0[];
    foreach() {
        if (ibm[] > 0 && ibm[] < 1 && (cr[] <= 1e-6 || cr[] >= ibm[] - 1e-6) &&
            !is_interior_cell(point, ibm, cr) && level == depth())
        {
            int near = 0;
            foreach_neighbor() //TODO: this is probably not necessary
                if (c[] > 0 && ibm[]) { near = 1; break; }
            ghostInter[] = near;
        }
        else
            ghostInter[] = 0;
    }

    foreach() {

        inter[] = c[] > 0 && c[] < ibm[];
        extra[] = inter[] && ibm[] > 0 && ibm[] < 1 && contact_angle[];

#if 1 // Make sure that extrapolation cells isn't an interior cell
        if (extra[]) {
            bool caCell = false;
            foreach_neighbor() {
                if (ibm[] && cr[]/ibm[] < 0.9) {
                  caCell = true;
                  break;
                }
            }
            extra[] = caCell;
        }
#endif

        if (extra[]) {
            coord ns = facet_normal (point, ibm, ibmf);
            double alphas = line_alpha (ibm[], ns);
            coord nf;

            foreach_dimension()
                nf.x = n.x[];

            normalize (&ns);
            normalize (&nf);
            coord nc = normal_contact (ns, nf, contact_angle[]);
            normalize_sum(&nc);

            foreach_dimension() 
                n.x[] = nc.x;

            alpha[] = immersed_alpha(c[], ibm[], nc, alpha[], ns, alphas, cr[]);
            if (extra[]) {
                c[] = plane_volume (nc, alpha[]);
                if (c[] <= 0) { // so ch and cr are consistent! prevents sudden removal of gginter
                    cr[] = 0;
                    if (!is_interior_cell(point, ibm, cr) && level == depth()) {
                        int near = 0;
                        foreach_neighbor() 
                            if (c[] > 0 && ibm[]) { near = 1; break; }
                        ghostInter[] = near;
                    }
                }
                inter[] = c[] > 0 && c[] < 1;
                extra[] = inter[] && ibm[] > 0 && ibm[] < 1 && !ghostInter[] && cr[] < ibm[]-1e-6 && contact_angle[];
            }
     }
    }

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
void set_contact_angle (scalar c, scalar cr0, const scalar ibm,
                        vector nf, scalar alphaf, vector ns, scalar alphas) // clean up names, more consistent (c -> ch)
{
    scalar ghostInter[];

    foreach() { gf00[] = c[]; }

    reconstruction_contact(c, cr0, nf, alphaf, ns, alphas, inter, ghostInter, extra);
    scalar c0[], c1[];
    scalar_clone (c0, c);
    scalar_clone (c1, c);

    #if 0
    foreach() {
        gginter[] = ghostInter[];
        c0[] = c[];
        gf0[] = c0[];
    }
    boundary({c0});
    #endif

    // 3. look for ghost and pure solid cells near triple point to enforce B.C
    foreach() {
        c0[] = c[];
        if (ibm[] <= 0. && level == depth()) {  // pure solid cell
            double ghostf = 0., totalWeight = 0.;
            coord ghostCell = {x,y,z};

            // 3.a. Calculate weights
            foreach_neighbor() {
                if (extra[]) {
                    // why fr0? to mitigate problems for a very slow moving contact line and
                    // a new cell only gets a little bit of f (due to the low velocity/flux)
                    double cellWeight = ibm[] * (1. - ibm[]) * cr0[] * (ibm[] - cr0[]);

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
                    if (ibm[] && cr0[] && !extra[])
                        check = true;
                    else if (ibm[] && extra[] ) {
                        check = false;
                        break;
                    }
                }
                if (check)
                    c[] = 1;
                else {
                    foreach_neighbor() {
                        if (ibm[] && cr0[] && !extra[])
                            check = true;
                        else if (ibm[] && extra[]) {
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

            // 3.c. extrapolate f from direct neighbors to get initial guess
            int cond0 = cr0[] <= 1e-6;
            coord ns0 = {-ns.x[], -ns.y[], -ns.z[]}, nf1;

            foreach_neighbor() {
                if (extra[]) {
                    nf1 = (coord){nf.x[], nf.y[], nf.z[]};

                    if ((contact_angle[] > 3.*pi/4. || contact_angle[] < pi/4.) &&
                        cond0 && dot_product_norm(nf1, ns0) <= 0)
                        continue;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }
                    
                    double cellWeight = ibm[] * (1. - ibm[]) * cr0[] * (ibm[] - cr0[]);

                    totalWeight += cellWeight;

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                }
            }

            if (totalWeight > 0.) {
                c[] = ghostf / totalWeight;
                if (c[] > 1 - INT_TOL) c[] = 1; // is this good?
            }
            // 3.e otherwise, cell should not be used to enforce C.A.
            else if (ghostf <= 0 && cr0[] <= 1e-6 ) {
                c[] = 0;
            }
        }
        c1[] = c[];
    }

#if 0
    scalar ctmp[];
    foreach() 
        ctmp[] = c[];
    boundary({ctmp,cr0,c});
#endif

    /**
    TODO: is this necessary for ALL contact angles? or just super hydrophobic?*/
    foreach() {
        gf1[] = c[];
        // make full cells with fractional ch conserve cr
        if (ghostInter[] && cr0[] >= ibm[] - 1e-6 && c[] != c0[] && c[] != 1.) {
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
                                         (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5})*ibm[];
            
            if (!approx_equal_double(fr, cr0[])) {
                const tripoint tcell = fill_tripoint (cr0[], mf, alpha, ns1, alphas[], c[], ibm[]);
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
            double alpha = immersed_alpha(c[], ibm[], nf1, alphaf[], ns1, alphas[], cr0[]);
            c[] = plane_volume(nf1, alpha);
        }
        #endif
    }

    clean_fluid (c, ibm);
    boundary ({c, cr0});
}
