
#define PRINTCA 1

extern scalar f;

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
    if (- ns.x * nf.y + ns.y * nf.x > 0) { // 2D cross product
        n.x = - ns.x * cos(angle) + ns.y * sin(angle);
        n.y = - ns.x * sin(angle) - ns.y * cos(angle);
    }
    else {
        n.x = - ns.x * cos(angle) - ns.y * sin(angle);
        n.y =   ns.x * sin(angle) - ns.y * cos(angle);
    }
    return n;
}


/*
This function extends the normal reconstruction() function in fractions.h by
taking into account the contact angle on immersed boundaries. Given a volume fraction
field, f, it fills fields n and alpha with the each cells normal and intercept, resp.
*/

void reconstruction_contact (scalar f, scalar fr, vector n, scalar alpha, scalar inter)
{
    reconstruction (f, n, alpha);

    foreach() {
        if (on_interface(ibm) && on_interface(f) && fr[]) {
            coord ns = facet_normal (point, ibm, ibmf);
            normalize (&ns);
            coord nf;

            foreach_dimension()
                nf.x = n.x[];
            coord nc = normal_contact (ns, nf, contact_angle[]);
            double mag = fabs(nc.x) + fabs(nc.y) + SEPS;
            foreach_dimension() {
                nc.x /= mag;
                n.x[] = nc.x;
            }
            alpha[] = line_alpha (f[], nc);
        }
        #if 0
        inter[] = real_interfacial (point, fr, f) && 
            !is_interior_cell (point, ibm, f, fr)? 1: 0;
        #else
        inter[] = on_interface(ibm) && f[] && 
                  !is_interior_cell (point, ibm, f, fr)? 1: 0;
        #endif
    }
    boundary ({n, alpha, inter});
}


/*
clean_fluid is used to remove any fluid inside the solid boundary that is
not necessary in enforcing the contact angle.
*/

void clean_fluid (scalar f, scalar fr, scalar ibm)
{
    foreach() {
        if (ibm[] == 0. && f[]) {
            int fluidNeighbors = 0;
            foreach_neighbor() {
                if (f[] > 0 && ibm[]) {
                    ++fluidNeighbors;
                    break;
                }
            }
            if (!fluidNeighbors) 
                f[] = 0;
        }
        #if 0
        if (on_interface(ibm) && f[] && !fr[]) {
            int realNeighbors = 0;
            foreach_neighbor() {
                if (fr[] > 0) {
                    ++realNeighbors;
                    break;
                }
            }
            if (!realNeighbors)
                f[] = 0;
        }
        #endif
    }
}


/*
This is the event that modifies the ghost cell's volume fraction, ghostf, at each timestep
necessary to impose the desired contact angle while conserving the real fluid, fr.

The event is named "tracer_advection" just to make sure that it is ran before the
surface tension force calculation which takes place in the acceleration event.
*/


double cerror3 = 0;

scalar inter[];

void set_contact_angle (scalar f, const scalar fr0, scalar ibm)
{
    #if PRINTCA
    fprintf (stderr, "\n=== SETTING C.A. ===\n");
    #endif

    scalar alphaf[], alphas[];
    vector nf0[], nf[], ns[];
    
    reconstruction (f, nf0, alphaf);
    reconstruction (ibm, ns, alphas);

    // 1. get real portion of fluid, fr (now a func. parameter)

    // 2. reconstruct n and alpha considering the contact angle B.C.
    reconstruction_contact (f, fr0, nf, alphaf, inter);

    // 3. look for ghost and pure solid cells near triple point to enforce B.C
    foreach() {
        if (ibm[] <= 0.) {  // pure solid cell
            double ghostf = 0., totalWeight = 0., alphaf1; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z}, nf2;
            int count = 0;

            // 3.a. Calculate weights
            foreach_neighbor() {

                if (on_interface(ibm) && inter[] && on_interface(f)) {
                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]); // changed from fr[]
                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf1 = {nf.x[], nf.y[]};
                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);

                    nf2 = nf1;
                    alphaf1 = alphaf[];
                    count++;
                }
            }
            
            // 3.b. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.) {
                #if PRINTCA && 1
                fprintf (stderr, "0th %d (%g, %g) nf ={%g,%g} f[]=%0.15g fr0[]=%0.15g tw=%g gf=%g ", 
                                  count, x, y, nf2.x, nf2.y, f[], fr0[], totalWeight, ghostf);
                fprintf (stderr, "gf/tw=%g alphaf1=%g\n", ghostf/totalWeight, alphaf1);
                #endif
                f[] = ghostf / totalWeight;
            }
        }
        #if 1
        else if (on_interface(ibm) && 
            (fr0[] <= INT_TOL || fr0[] >= ibm[] - INT_TOL)) // solid interface cell with full
        {                                                                           // or empty real fluid
            double ghostf = 0., totalWeight = 0., alphaf1 = 0;
            coord ghostCell = {x,y,z}, nf2;
            int count = 0;

            // 3.c. extrapolate f from direct neighbors to get initial guess
            foreach_neighbor() {
                if (on_interface(ibm) && on_interface(f) && inter[]) {

                    double cellWeight = ibm[] * (1. - ibm[]) * fr0[] * (1. - fr0[]);
                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf1 = {nf.x[], nf.y[]};
                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);

                    nf2 = nf1;
                    alphaf1 = alphaf[];
                    count++;
                }
            }

            // 3.d. if extrapolation results in f inside the interface cell, adjust it to preserve fr
            if (totalWeight > 0. && ghostf > 0.) {
                coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};

                #if PRINTCA
                fprintf (stderr, "1st %d (%g, %g) nf ={%g,%g} f[]=%0.15g fr0[]=%0.15g tw=%g gf=%g ", 
                                  count, x, y, nf2.x, nf2.y, f[], fr0[], totalWeight, ghostf);
                fprintf (stderr, "gf/tw=%g alphaf1=%g\n", ghostf/totalWeight, alphaf1);
                #endif

                f[] = ghostf / totalWeight;
                if (f[] > 1 - INT_TOL) f[] = 1;

                #if PRINTCA
                fprintf (stderr, "  f`[]=%0.15g \n", f[]);
                #endif
            }
            // 3.e otherwise, cell should not be used to enforce C.A.
            else if (ghostf <= 0 && fr0[] <= 0)
                f[] = 0;
        }
        #endif
    }

    // 4. Correct f in cells that may violate volume conservation by calculating
    //    and comparing fr before and after step 3.

    #if 1
    boundary ({f});
    reconstruction (f, nf, alphaf);

    scalar fr[];
    real_fluid (f, fr); // should be combined with following foreach loop for efficiency

    foreach() {
          if (fr[] < fr0[] - 1e-10 || fr[] > fr0[] + 1e-10) {

            if (fr0[] >= 1-1e-10) {
                f[] = 1;
                continue;
            }

            if (!ns.x[] && !ns.y[]) {
                f[] = fr0[];
                continue;
            }

            coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
            double freal = ibm[]*immersed_fraction (f[], nf1, alphaf[], ns1, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);

            if (fabs(freal - fr0[]) <= 1e-6)
                continue;

            if (!nf1.x && !nf1.y) {
                fprintf (stderr, "WARNING: liquid normal vector is zero in the CA function!\n");
                f[] = 0.99999;
                nf1 = interface_normal (point, f);
            }

            #if PRINTCA
            fprintf (stderr, "2nd (%g, %g) nf={%g,%g} ns={%g,%g} f[]=%0.15g fr0[]=%0.15g fr[]=%0.15g freal=%g ibm=%g\n", 
                             x, y, nf1.x, nf1.y, ns1.x, ns1.y, f[], fr0[], fr[], freal, ibm[]);
            #endif

            coord lhs = {-0.5,-0.5}, rhs = {0.5,0.5}, ip[2]; // intersecting point
            int bp = boundary_points (ns1, alphas[], lhs, rhs, ip);
            double alpha0 = nf1.x*ip[0].x + nf1.y*ip[0].y;
            double alpha1 = nf1.x*ip[1].x + nf1.y*ip[1].y;

            double f0 = plane_volume (nf1, alpha0);
            double f1 = plane_volume (nf1, alpha1);

            double freal0 = ibm[]*immersed_fraction (f0, nf1, alpha0, ns1, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5});
            double freal1 = ibm[]*immersed_fraction (f1, nf1, alpha1, ns1, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5});

            #if PRINTCA
            fprintf(stderr, "    f0=%g f1=%g freal0=%g freal1=%g alpha0=%g alpha1=%g\n",
                f0, f1, freal0, freal1, alpha0, alpha1);
            fprintf(stderr, "    ip0 = {%g, %g} ip1 = {%g, %g}\n",
                ip[0].x, ip[0].y, ip[1].x, ip[1].y);
            #endif

            double alpha = alphaf[];
            if (fabs(freal0 - fr0[]) <= 1e-8)
                alpha = alpha0;
            else if (fabs(freal1 - fr0[]) <= 1e-8)
                alpha = alpha1;
            else {
                fr0[] = clamp(fr0[], 0., ibm[]);
                tripoint tcell = fill_tripoint (fr0[], nf1, alphaf[], ns1, alphas[], f[], ibm[]);
                bisection_solver (&alpha, -0.75, 0.75, tcell);
            }

            f[] = plane_volume (nf1, alpha);
            if (f[] > 1 - INT_TOL) f[] = 1;
            #if PRINTCA
            fprintf (stderr, "  f`[]=%0.15g \n", f[]);
            #endif
        }
    }
    #endif

   cerror3 = real_volume(f);

    clean_fluid(f, fr0, ibm);
    boundary ({f});
}
