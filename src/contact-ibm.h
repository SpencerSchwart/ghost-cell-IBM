#define CA 1

#include "ibm-gcm-vof.h"

#undef INT_TOL
#define INT_TOL 1e-10

/*
normal_contact returns the properly-oriented normal of an interface touching the 
immersed boundary. ns is the normal to the immersed solid boundary and nf is the
normal to the VOF interface not taking into account the contact angle boundary
condition. angle is the imposed static contact angle.

TODO: show derivation better
TODO: 3D extension
*/

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
#endif

    return n;
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
field, f, it fills fields n and alpha with the each cells normal and intercept, resp. */

void reconstruction_contact (scalar c, scalar cr, vector n, scalar alpha, 
                                  vector ns, scalar alphas, scalar inter, 
                                  scalar ghostInter, scalar extra)
{
    scalar c0[];
    foreach() {
        if (ibm[] > 0 && ibm[] < 1 && (cr[] <= 1e-10 || cr[] >= ibm[] - 1e-10) &&
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

        inter[] = c[] > 0 && c[] < 1;
        extra[] = inter[] && ibm[] > 0 && ibm[] < 1;

#if 0
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
            coord nc = normal_contact (ns, nf, contact_angle[]);

            double mag = 0;
            foreach_dimension()
                mag += fabs(nc.x);

            foreach_dimension() {
                nc.x /= mag + SEPS;
                n.x[] = nc.x;
            }
            alpha[] = immersed_alpha(c[], ibm[], nc, alpha[], ns, alphas, cr[]);
            if (extra[]) {
                c[] = plane_volume (nc, alpha[]);
                if (c[] <= 0) { // so ch and cr are consistent! prevents sudden removal of gginter
                    cr[] = 0;
                    if (!is_interior_cell(point, ibm, cr) && level == depth()) {
                        int near = 0;
                        foreach_neighbor() //TODO: this is probably not necessary
                            if (c[] > 0 && ibm[]) { near = 1; break; }
                        ghostInter[] = near;
                    }
                }
                inter[] = c[] > 0 && c[] < 1;
                extra[] = inter[] && ibm[] > 0 && ibm[] < 1 && !ghostInter[] && cr[] < ibm[]-1e-10;
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
                        vector nf, scalar alphaf, vector ns, scalar alphas)
{
    scalar ghostInter[];

    //foreach() { gf00[] = c[]; }

    reconstruction_contact(c, cr0, nf, alphaf, ns, alphas, inter, ghostInter, extra);
    scalar c0[];
    foreach() {
        gginter[] = ghostInter[];
        c0[] = c[];
        gf0[] = c0[];
    }
    boundary({c0});

    // 3. look for ghost and pure solid cells near triple point to enforce B.C
    foreach() {
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
            else if (contact_angle[] <= 0.5*pi) {
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
        }
        if (ghostInter[])
        {
            double ghostf = 0., totalWeight = 0.;
            coord ghostCell = {x,y,z};

            // 3.c. extrapolate f from direct neighbors to get initial guess
            int cond0 = cr0[] <= 1e-10;
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
            else if (ghostf <= 0 && cr0[] <= 1e-10 ) {
                c[] = 0;
            }
        }
    }

    scalar ctmp[];
    foreach() 
        ctmp[] = c[];
    boundary({ctmp,cr0,c});

    foreach() {
        gf1[] = c[];
        // make full cells with fractional ch conserve cr
        if (ghostInter[] && cr0[] >= ibm[] - 1e-10 && c[] != c0[] && c[] != 1.) {
            coord mf = interface_normal(point, ctmp);
            if (!mf.x && !mf.y && !mf.z) {
                double ctemp = ctmp[];
                ctmp[] = 0.5;
                mf = interface_normal(point, ctmp);
                ctmp[] = ctemp;
            }
            coord ns1 = {ns.x[], ns.y[], ns.z[]};
            double alpha = plane_alpha(ctmp[], mf);           
            double fr = immersed_fraction(ctmp[], mf, alpha, ns1, alphas[], 
                                         (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5})*ibm[];
            
            if (!approx_equal_double(fr, cr0[])) {
                const tripoint tcell = fill_tripoint (cr0[], mf, alpha, ns1, alphas[], c[], ibm[]);
                double alphatmp = ghost_alpha(tcell, -0.6, 0.6);
                c[] = plane_volume(mf, alphatmp);
            }
        }
        if (extra[]) {
            coord nf1 = interface_normal(point, ctmp);
            if (!nf1.x && !nf1.y && !nf1.z) {
              c[] = 0;
              continue;
            }
            coord ns1 = {ns.x[], ns.y[], ns.z[]};
            double alpha = immersed_alpha(c[], ibm[], nf1, alphaf[], ns1, alphas[], cr0[]);
            c[] = plane_volume(nf1, alpha);
        }
    }

    clean_fluid (c, ibm);
    boundary ({c});
}
