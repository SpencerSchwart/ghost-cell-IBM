
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
This function returns true if there exists a triple points within the cell, 
i.e. the liquid interface intersects the solid one, given the normals of the
liquid and solid interfaces.
*/

bool is_triple_point (Point point, coord nf, coord ns)
{
    if (!(on_interface(ibm)) || !(on_interface(f)))
        return false;
    if (ns.x == 0 || ns.y == 0 || nf.x == 0 || nf.y == 0)
        return false;

    double alphas = plane_alpha (ibm[], ns);

    double alphaf = plane_alpha (f[], nf);

    double intercept = ((alphas/ns.y) - (alphaf/nf.y)) /
                       ((ns.x/ns.y) - (nf.x/nf.y));
    
    return fabs(intercept) <= 0.5;

}


/*
This function extends the normal reconstruction() function in fractions.h by
taking into account the contact angle on immersed boundaries. Given a volume fraction
field, f, it fills fields n and alpha with the each cells normal and intercept, resp.
*/

void reconstruction_contact (scalar f, scalar fr, vector n, scalar alpha)
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
    }
    boundary ({n, alpha});
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
    }
}


/*
This is the event that modifies the ghost cell's volume fraction, ghostf, at each timestep
necessary to impose the desired contact angle while conserving the real fluid, fr.

The event is named "tracer_advection" just to make sure that it is ran before the
surface tension force calculation which takes place in the acceleration event.
*/

event tracer_advection (i++)
{
    scalar alphaf[], alphas[];
    vector nf0[], nf[], ns[];
    
    reconstruction (f, nf0, alphaf);
    reconstruction (ibm, ns, alphas);

    // 1. get real portion of fluid, fr
    scalar fr0[];
    real_fluid (f, fr0);

    // 2. reconstruct n and alpha considering the contact angle B.C.
    reconstruction_contact (f, fr0, nf, alphaf);

    // 3. look for ghost and pure solid cells near triple point to enforce B.C
    foreach() {
        if (ibm[] <= 0.) {  // pure solid cell
            double ghostf = 0., totalWeight = 0.; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z};
    
            // 3.a. Calculate weights
            foreach_neighbor() {
                if (on_interface(ibm) && (on_interface(fr0) &&
                    fr0[] <= ibm[] - INT_TOL && fr0[] > 0 + INT_TOL)) { // cell with potential triple point

                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);
                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf1 = {nf.x[], nf.y[]};
                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                }
            }
            
            // 3.b. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.)
                f[] = ghostf / totalWeight;
        }
        else if (on_interface(ibm) && (fr0[] <= 0+INT_TOL || fr0[] >= ibm[]-INT_TOL)) // solid interface cell with full
        {                                                           // or empty real fluid
            double ghostf = 0., totalWeight = 0., alphaf1 = 0;
            coord ghostCell = {x,y,z}, nf2;
            int count = 0;

            // 3.c. extrapolate f from direct neighbors to get initial guess
            foreach_neighbor(1) {
                if (on_interface(ibm) && on_interface(f) && fr0[] < ibm[] && fr0[] > 0 + INT_TOL) {

                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);
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
#if 1
                coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
                fprintf (stderr, "1st %d (%g, %g) nf ={%g,%g} f[]=%0.15g fr0[]=%0.15g tw=%g gf=%g ", 
                                  count, x, y, nf2.x, nf2.y, f[], fr0[], totalWeight, ghostf);
                fprintf (stderr, "gf/tw=%g alphaf1=%g\n", ghostf/totalWeight, alphaf1);
                //if (ghostf / totalWeight >= 1. && fr0[] >= ibm[] - INT_TOL) {
                if (ghostf / totalWeight >= 1.) {
                    fprintf(stderr, "(a) fr0[]=%g ibm[]=%g\n", fr0[], ibm[]);     
                    f[] = ghostf / totalWeight;
                }
                else { // should it really be nf2 here instead of nf1?
                    fprintf(stderr, "(b) fr0[]=%g ibm[]=%g tw=%g gf=%g\n", fr0[], ibm[], totalWeight, ghostf);
                    //double alpha = immersed_line_alpha (point, nf2, alphaf1, ns1, alphas[], fr0[]);
                    //f[] = plane_volume (nf2, alpha);
                    f[] = ghostf / totalWeight;
                }
                if (f[] > 1 - INT_TOL) f[] = 1;
#endif
                fprintf (stderr, "  f`[]=%0.15g \n", f[]);
            }
            // 3.e otherwise, cell should not be used to enforce C.A.
            else if (ghostf <= 0 && fr0[] <= 0)
                f[] = 0;
        }
    }

    // 4. Correct f in cells that may violate volume conservation by calculating
    //    and comparing fr before and after step 3.

    boundary ({f});
    reconstruction (f, nf, alphaf);

    scalar fr[];
    real_fluid (f, fr);

    foreach() {
          if (fr[] < fr0[] - INT_TOL || fr[] > fr0[] + INT_TOL) {

            coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
            double freal = ibm[]*immersed_fraction (f[], nf1, alphaf[], ns1, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);

            if (freal >= fr0[] - 1e-10 && freal <= fr0[] + 1e-10)
                continue;

            fprintf (stderr, "2nd (%g, %g) nf ={%g,%g} f[]=%0.15g fr0[]=%0.15g fr[]=%0.15g freal=%g\n", 
                             x, y, nf1.x, nf1.y, f[], fr0[], fr[], freal);
            double alpha = immersed_line_alpha (point, nf1, alphaf[], ns1, alphas[], fr0[]);
            f[] = plane_volume (nf1, alpha);
            if (f[] > 1 - INT_TOL) f[] = 1;
            fprintf (stderr, "  f`[]=%0.15g \n", f[]);
        }
    }
    clean_fluid(f, fr0, ibm);
    boundary ({f});
}
