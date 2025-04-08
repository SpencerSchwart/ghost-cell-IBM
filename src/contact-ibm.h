
//(const) scalar contact_angle;
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
#endif

/*
get_normal_contact returns the normal of the fluid considering the contact angle
boundary condition and the orientation of the solid interface.

Note: both ns and nf must be normalized so |n.x| + |n.y| = 1,
      the returned coord, n, is also normalized in such a way
*/
#if 0
static inline coord normal_contact (coord ns, coord nf, double angle)
{
    double norm = distance(ns.x,ns.y) + SEPS;
    ns.x /= norm, ns.y /= norm;
    norm = distance(nf.x,nf.y);
    nf.x /= norm, nf.y /= norm;

    coord n = {0,0};
    if (- ns.x * nf.y + ns.y * nf.x > 0) { // 2D cross product
        n.x = - ns.x * cos(angle) + ns.y * sin(angle);
        n.y = - ns.x * sin(angle) - ns.y * cos(angle);
    }
    else {
        n.x = - ns.x * cos(angle) - ns.y * sin(angle);
        n.y =   ns.x * sin(angle) - ns.y * cos(angle);
    }

    norm = abs(n.x) + abs(n.y) + SEPS;
    n.x /= norm, n.y /= norm;

    return n;
}
#endif

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
            //fprintf (stderr, "| (%g, %g) nc={%g, %g} alpha=%g\n", x, y, nc.x, nc.y, alpha[]);
        }
    }

    boundary ({n, alpha});
}

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
            if (!fluidNeighbors) {
                f[] = 0;
            }
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
necessary to impose the desired contact angle.

TODO: make this event run before surface tension (acceleration) event
*/
#if 1
scalar fr[];
event tracer_advection (i++)
{
    scalar alphaf[], alphas[];
    vector nf0[], nf[], ns[];
    
    reconstruction (f, nf0, alphaf);
    reconstruction (ibm, ns, alphas);

    // 2. get real portion of fluid, fr
    real_fluid (f, fr);

    // 1. reconstruct n and alpha considering the contact angle B.C.
    reconstruction_contact (f, fr, nf, alphaf);

    // 3. look for ghost or pure solid cells near triple point to enforce B.C
    foreach() {
        if (ibm[] <= 0.) {
            double ghostf = 0., totalWeight = 0.; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z};
    
            // 3. Calculate weights
            foreach_neighbor() {
              if (on_interface(ibm) && (on_interface(f) && fr[] <= ibm[] - 1e-6 && fr[] > 1e-6)) { // cell with potential triple point
                coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf0.x[], nf0.y[]};
              //if (is_triple_point(point, nf1, ns1)) {
                coord nftemp, nstemp;
                foreach_dimension() {
                    nftemp.x = nf.x[];
                    nstemp.x = ns.x[];
                }
                double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);

                totalWeight += cellWeight;

                coord leftPoint = {x, y, z}, rightPoint;
                foreach_dimension() {
                    leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                    rightPoint.x = leftPoint.x + 1.;
                }

                coord nf2;
                foreach_dimension()
                    nf2.x = nf.x[];

                ghostf += cellWeight * rectangle_fraction (nf2, alphaf[], leftPoint, rightPoint);
                    
                }
            }
            
            // 4. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.) {
                f[] = ghostf / totalWeight;
                //fprintf (stderr, "ghost fluid (%g, %g) = %g, gf=%g tw=%g\n", x, y, f[], ghostf, totalWeight);
            }
        }
        #if 1
        #if 1
        else if (on_interface(ibm) && (fr[] <= 0 || fr[] >= ibm[]))
        {
            double ghostf = 0., totalWeight = 0., alphaf1 = 0;
            coord ghostCell = {x,y,z}, nf2;

            int count = 0;
            foreach_neighbor(1) {
                if (on_interface(ibm) && on_interface(f) && fr[] < ibm[] && fr[] > 0) { // fr > 0? or fr < ibm?
                    coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);

                    totalWeight += cellWeight;
                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                    nf2 = nf1;
                    alphaf1 = alphaf[];
                    count++;
                }
            }
            //fprintf (stderr, "(%g, %g) fr=%g tw=%0.15g gf=%0.15g\n", x, y, fr[], totalWeight, ghostf);
            if (totalWeight > 0. && ghostf > 0. && f[] >= 0) {
                
                coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
                //fprintf (stderr, "1st %d (%g, %g) nf ={%g,%g} f[]=%0.15g fr[]=%0.15g tw=%g gf=%g ", 
                //                  count, x, y, nf2.x, nf2.y, f[], fr[], totalWeight, ghostf);
                //fprintf (stderr, "gf/tw=%g alphaf1=%g\n", ghostf/totalWeight, alphaf1);
                if (ghostf / totalWeight >= 1.)
                     f[] = ghostf / totalWeight;
                else {
                    double alpha = immersed_line_alpha (point, nf2, alphaf1, ns1, alphas[], fr[]);
                    f[] = plane_volume (nf2, alpha);
                }
                if (f[] > 1-1e-6)
                    f[] = 1;
                //fprintf (stderr, "  f`[]=%0.15g\n", f[]);
            }
            else if (ghostf <= 0 && fr[] <= 0) {
                //fprintf(stderr, "empty cell\n");
                f[] = 0;
            }
        }
        #else // works for incline
        else if (on_interface(ibm) && fr[] <= 0)
        {
            double ghostf = 0., totalWeight = 0.;
            coord ghostCell = {x,y,z}, nf2;

            int count = 0;
            foreach_neighbor() {
                if (on_interface(ibm) && on_interface(f)) { // fr > 0?
                    coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);

                    totalWeight += cellWeight;
                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                    nf2 = nf1;
                    count++;
                }
            }

            if (totalWeight > 0. && f[] >= 0) {
                coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
                //fprintf (stderr, "%d (%g, %g) nf ={%g,%g} f[]=%g fr[]=%g\n", 
                //                  count, x, y, nf2.x, nf2.y, f[], fr[]);
                double alpha = immersed_line_alpha (point, nf2, alphaf[], ns1, alphas[], fr[]);
                //f[] = !plane_volume (nf1, alpha)? ghostf/ totalWeight: plane_volume (nf1,alpha);
                f[] = ghostf / totalWeight;
                //fprintf (stderr, "  f`[]=%g\n", f[]);
            }
        }
        #endif
        #endif
    }

    // Correct f in cells that may violate volume conservation
#if 1
    reconstruction (f, nf, alphaf);

    foreach() {
        if (on_interface(ibm) && on_interface(f) && (fr[] <= 0 || fr[] >= ibm[])) {
               
            coord ns1 = {ns.x[], ns.y[]}, nf1 = {nf.x[], nf.y[]};
            double freal = ibm[]*immersed_fraction (f[], nf1, alphaf[], ns1, alphas[],
                                              (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5},0);
            if (freal == fr[])
                continue;

            //fprintf (stderr, "2nd (%g, %g) nf ={%g,%g} f[]=%g fr[]=%g\n", 
            //                 x, y, nf1.x, nf1.y, f[], fr[]);
            double alpha = immersed_line_alpha (point, nf1, alphaf[], ns1, alphas[], fr[]);
            f[] = plane_volume (nf1, alpha);
            if (f[] > 1-1e-6)
                f[] = 1;
            //if (ghostf / totalWeight >= 1.)
            //     f[] = ghostf / totalWeight;
            //fprintf (stderr, "  f`[]=%g\n", f[]);
        }
    }
#endif
    clean_fluid(f, fr, ibm);
    boundary ({f});
}
#endif



#if 0
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
    //clean_fluid(f, ibm);
    boundary ({f});
}
#endif
