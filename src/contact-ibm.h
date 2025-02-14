
(const) scalar contact_angle;


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

void reconstruction_contact (scalar f, vector n, scalar alpha)
{
    reconstruction (f, n, alpha);

    foreach() {
        if (on_interface(ibm) && on_interface(f)) {
            coord ns = facet_normal (point, ibm, ibmf);
            normalize (&ns);
            coord nf;

            foreach_dimension()
                nf.x = n.x[];
            coord nc = normal_contact (ns, nf, contact_angle[]);

            foreach_dimension()
                n.x[] = nc.x;
            alpha[] = line_alpha (f[], nc);
        }
    }

    boundary ({n, alpha});
}


/*
clean_fluid is used to remove any fluid inside the solid boundary that are
not necessary in enforcing the contact angle.
*/

void clean_fluid (scalar f, scalar ibm)
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
            if (!fluidNeighbors) {
                f[] = 0;
            }
        }
    }
}

/*
This is the event that modifies the ghost cell's volume fraction, ghostf, at each timestep
necessary to impose the desired contact angle.
*/

event contact (i++)
{
    vector n[];
    scalar alpha[];
    
    // 1. reconstruct n and alpha considering the contact angle B.C.
    reconstruction_contact (f, n, alpha);

    // 2. look for ghost or pure solid cells near triple point to enforce B.C
    foreach() {
        if (ibm[] == 0.) {
            double ghostf = 0., totalWeight = 0.; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z};
    
            // 3. Calculate weights
            foreach_neighbor() {
                if (on_interface(ibm) && on_interface(f)) { // cell with potential triple point
                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);
                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf;
                    foreach_dimension()
                        nf.x = n.x[];

                    ghostf += cellWeight * rectangle_fraction (nf, alpha[], leftPoint, rightPoint);
                    
#if 0
                    fprintf (stderr, "cw=%g tw=%g nf.x=%g nf.y=%g ghostf=%g ibm=%g f=%g\n",
                                      cellWeight, totalWeight, nf.x, nf.y, ghostf, ibm[], f[]);
                    fprintf (stderr, "lp.x=%g lp.y=%g rp.x=%g rp.y=%g alpha=%g\n",
                                      leftPoint.x, leftPoint.y, rightPoint.x, rightPoint.y, alpha[]);
#endif                                      
                }
            }
            
            // 4. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.) {
                f[] = ghostf / totalWeight;
#if 0                
                fprintf(stderr, "f[]=%g ghostf=%g tw=%g\n", f[], ghostf, totalWeight);
                fprintf (stderr, "### Solid Cell @(%g,%g)###\n\n", x, y);
#endif                
            }
        }
    }
    clean_fluid(f, ibm);
    boundary ({f});
}

