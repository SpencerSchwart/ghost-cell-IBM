
(const) scalar contact_angle;
//extern scalar contact_angle;
extern scalar f;

/*
normal_contact returns the properly-oriented normal of an interface touching the 
immersed boundary. ns is the normal to the immersed solid boundary and nf is the
normal to the VOF interface not taking into account the contact angle boundary
condition. angle is the imposed static contact angle.

TODO: show derivation better
TODO: 3D extension
*/

#if 1
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


#if 1
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
#endif

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
#endif

void reconstruction_contact_vof (scalar f, vector n, scalar alpha)
{
    reconstruction (f, n, alpha);

    vector normals[];
    foreach() {
        if (on_interface(ibm)) {
            coord ns, midPoint;
            double alphaSolid = 0;
            ibm_geometry (point, &ns, &midPoint, &alphaSolid);
            double mag = fabs(ns.x) + fabs(ns.y);
            ns.x /= mag, ns.y /= mag;
            foreach_dimension()
                normals.x[] = ns.x;
        }

    }

    foreach() {
       if (on_interface(ibm) && on_interface(f)) {
        coord ns1 = {normals.x[], normals.y[]}, nf1 = {n.x[], n.y[]};
        //if (is_triple_point(point, ns1, nf1)) {
            coord ns = facet_normal (point, ibm, ibmf);
            normalize (&ns);
            coord nf;

            foreach_dimension()
                nf.x = n.x[];
            coord nc = normal_contact (ns, nf, contact_angle[]);
            //double mag = fabs(nc.x) + fabs(nc.y);
            double mag = 1;
            foreach_dimension() 
                n.x[] = nc.x/mag;

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
#if 1
event contact (i++)
{
    vector n[];
    scalar alpha[];

    vector nf0[];
    
    reconstruction (f, nf0, alpha);

    // 1. reconstruct n and alpha considering the contact angle B.C.
    reconstruction_contact (f, n, alpha);

    vector normals[];
    scalar alphas[];
    foreach() {
        if (on_interface(ibm)) {
            coord ns, midPoint;
            double alphaSolid = 0;
            ibm_geometry (point, &ns, &midPoint, &alphaSolid);
            double mag = fabs(ns.x) + fabs(ns.y);
            ns.x /= mag, ns.y /= mag;
            foreach_dimension()
                normals.x[] = ns.x;
            alphas[] = alphaSolid;
        }

    }

    // 2. look for ghost or pure solid cells near triple point to enforce B.C
    foreach() {
        if (ibm[] <= 0.) {
            double ghostf = 0., totalWeight = 0.; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z};
    
            // 3. Calculate weights
            foreach_neighbor() {
                if (on_interface(ibm) && on_interface(f)) { // cell with potential triple point
                coord ns = {normals.x[], normals.y[]}, nf1 = {nf0.x[], nf0.y[]};
                //if (is_triple_point(point, nf1, ns)) {
                coord nftemp, nstemp;
                foreach_dimension() {
                    nftemp.x = n.x[];
                    nstemp.x = normals.x[];
                }
                double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);
                //double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]) + is_triple_point(point,nf1, ns)/6;

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
                    fprintf (stderr, "lp.x=%g lp.y=%g rp.x=%g rp.y=%g alpha=%g x=%g y=%g\n",
                                      leftPoint.x, leftPoint.y, rightPoint.x, rightPoint.y, alpha[], x, y);
                    fprintf (stderr, "ns.x=%g ns.y=%g alphas=%g\n", normals.x[], normals.y[], alphas[]); 
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
#endif

void impose_contact_angle (scalar f, scalar ibm)
{
    vector n[];
    scalar alpha[];

    vector nf0[];
    reconstruction (f, nf0, alpha);

    reconstruction_contact (f, n, alpha);

    vector normals[];
    scalar alphas[];
    foreach() {
        if (on_interface(ibm)) {
            coord ns, midPoint;
            double alphaSolid = 0;
            ibm_geometry (point, &ns, &midPoint, &alphaSolid);
            double mag = fabs(ns.x) + fabs(ns.y);
            ns.x /= mag, ns.y /= mag;
            foreach_dimension()
                normals.x[] = ns.x;
        }
    }

    foreach() {
        if (ibm[] <= 0.5) {
            double ghostf = 0., totalWeight = 0.; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z};
    
            // 3. Calculate weights
            foreach_neighbor() {
                //if (on_interface(ibm) && on_interface(f)) { // cell with potential triple point
                coord ns = {normals.x[], normals.y[]}, nf1 = {nf0.x[], nf0.y[]};
                if (is_triple_point(point, nf1, ns)) {
                    coord nftemp, nstemp;
                    foreach_dimension() {
                        nftemp.x = n.x[];
                        nstemp.x = normals.x[];
                    }
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
                }
            }
            if (totalWeight > 0.)
                f[] = ghostf / totalWeight;
        }
    }
    clean_fluid(f, ibm);
    boundary ({f});
}
