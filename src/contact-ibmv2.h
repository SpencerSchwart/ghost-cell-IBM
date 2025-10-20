

#define PRINTCA 0

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


void clean_fluid (scalar f, scalar fr, scalar ibm)
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


/*
This function extends the normal reconstruction() function in fractions.h by
taking into account the contact angle on immersed boundaries. Given a volume fraction
field, f, it fills fields n and alpha with the each cells normal and intercept, resp.
*/


/*
clean_fluid is used to remove any fluid inside the solid boundary that is
not necessary in enforcing the contact angle.
*/

void reconstruction_contact_test (scalar f, scalar fr, vector n, scalar alpha, 
                                  vector ns, scalar alphas, scalar inter, 
                                  scalar ghostInter, scalar extra)
{
    scalar f0[];
    foreach() {
        if (ibm[] > 0 && ibm[] < 1 && (fr[] <= 1e-10 || fr[] >= ibm[] - 1e-10) &&
            !is_interior_cell(point, ibm, fr) && level == depth())
        {
            int near = 0;
            foreach_neighbor() //TODO: this is probably not necessary
                if (f[] > 0 && ibm[]) { near = 1; break; }
            ghostInter[] = near;
        }
        else
            ghostInter[] = 0;
    }

    foreach(serial) {

        inter[] = f[] < 1 && f[] > 0;
        extra[] = inter[] && ibm[] > 0 && ibm[] < 1 && !ghostInter[] && fr[] < ibm[]-1e-10;

        if (on_interface(ibm) && f[] > 1e-10 && f[] < 1 - 1e-10) {
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
            alpha[] = immersed_alpha_temp (f[], ibm[], nc, alpha[], ns, alphas, fr[]);
            if (extra[]) {
                f[] = plane_volume (nc, alpha[]);
                if (f[] <= 0) { // so ch and cr are consistent! prevents sudden removal of gginter
                    fr[] = 0;
                    if (!is_interior_cell(point, ibm, fr) && level == depth()) {
                        int near = 0;
                        foreach_neighbor() //TODO: this is probably not necessary
                            if (f[] > 0 && ibm[]) { near = 1; break; }
                        ghostInter[] = near;
                    }
                }
                inter[] = f[] < 1 && f[] > 0;
                extra[] = inter[] && ibm[] > 0 && ibm[] < 1 && !ghostInter[] && fr[] < ibm[]-1e-10;
            }
        }
    }

    boundary ({f, fr, n, alpha, inter, ghostInter, extra});
}

scalar gftest[];

void set_contact_angle_solid (scalar f, const scalar ibm, vector nf, scalar alphaf, 
                              vector ns, scalar alphas, scalar extra)
{
    foreach() {
        if (ibm[] <= 0. && level == depth()) {  // pure solid cell
            double ghostf = 0., totalWeight = 0., alphaf1; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z}, nf2;
            int count = 0;

            // 3.a. Calculate weights
            foreach_neighbor() {

                if (extra[]) {
                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);

                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf1 = {nf.x[], nf.y[], nf.z[]}, ns1 = {ns.x[], ns.y[], ns.z[]};

                    //double alphatmp = immersed_alpha_temp (f[], ibm[], nf1, alphaf[], ns1, alphas[], fr0[]);

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);

                    #if PRINTCA && 0
                    double gf = rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                    fprintf(stderr, " 1st pre (%g, %g, %g) f = %g ibm = %g gf = %g\n", 
                        x, y, z, f[], ibm[], gf);
                    #endif

                    nf2 = nf1;
                    alphaf1 = alphaf[];
                    count++;
                }
            }
            
            // 3.b. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.) 
                f[] = ghostf / totalWeight;
            else
                f[] = 0;
        }
        gftest[] = f[];
    }
}

/*
This is the event that modifies the ghost cell's volume fraction, ghostf, at each timestep
necessary to impose the desired contact angle while conserving the real fluid, fr.

The event is named "tracer_advection" just to make sure that it is ran before the
surface tension force calculation which takes place in the acceleration event.
*/


scalar gginter[];
scalar gf0[], gf1[], gf2[], gf3[], gf4[];

trace
void set_contact_angle_tension (scalar f, scalar fr0, const scalar ibm,
                                vector nf, scalar alphaf, vector ns, scalar alphas)
{
    #if PRINTCA
    fprintf (stderr, "\n=== SETTING TENSION C.A. ===\n");
    #endif
    scalar inter[];
    scalar ghostInter[];
    scalar extra[];

    // 1. get real portion of fluid, fr (now a func. parameter)

    // 2. reconstruct n and alpha considering the contact angle B.C.
    #if PRINTCA
    fprintf(stderr, "starting reconstruction contact test\n");
    #endif
    reconstruction_contact_test (f, fr0, nf, alphaf, ns, alphas, inter, ghostInter, extra);
#if 1
    scalar f0[];
    foreach() {
        gginter[] = ghostInter[];
        f0[] = f[];
        gf0[] = f0[];
    }
    boundary({f0});
#endif
    #if PRINTCA
    fprintf(stderr, "starting main contact loop\n");
    #endif

    // 3. look for ghost and pure solid cells near triple point to enforce B.C
    foreach() {
        #if 1
        if (ibm[] <= 0. && level == depth()) {  // pure solid cell
            double ghostf = 0., totalWeight = 0., alphaf1;
            coord ghostCell = {x,y,z}, nf2;
            int count = 0;

            // 3.a. Calculate weights
            foreach_neighbor() {
                if (extra[]) {
                    double cellWeight = ibm[] * (1. - ibm[]) * f0[] * (1. - f0[]);

                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf1 = {nf.x[], nf.y[], nf.z[]}, ns1 = {ns.x[], ns.y[], ns.z[]};

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);

                    #if PRINTCA && 0
                    double gf = rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);
                    fprintf(stderr, " 1st pre (%g, %g, %g) f = %g ibm = %g gf = %g\n", 
                        x, y, z, f[], ibm[], gf);
                    #endif

                    nf2 = nf1;
                    alphaf1 = alphaf[];
                    count++;
                }
            }
            
            // 3.b. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.) {
                #if PRINTCA && 0
                fprintf (stderr, "0th %d (%g, %g, %g) nf ={%g,%g} f[]=%0.15g fr0[]=%0.15g tw=%g gf=%g ", 
                                  count, x, y, z, nf2.x, nf2.y, f[], fr0[], totalWeight, ghostf);
                fprintf (stderr, "gf/tw=%g alphaf1=%g\n", ghostf/totalWeight, alphaf1);
                #endif
                f[] = ghostf / (totalWeight + SEPS);
            }
            #if 1
            else {
                bool check = false;
                foreach_neighbor(1) {
                    if (ibm[] && fr0[] >= ibm[] - 1e-6)
                        check = true;
                    else if (ibm[] && fr0[] < ibm[] - 1e-6) {
                        check = false;
                        break;
                    }
                }
                if (check) {
                    f[] = 1;
                }
                foreach_neighbor() {
                    if (ibm[] && fr0[] >= ibm[] - 1e-6)
                        check = true;
                    else if (ibm[] && fr0[] < ibm[] - 1e-6) {
                        check = false;
                        break;
                    }
                }
                if (check) {
                    f[] = 1;
                }
            }
            #endif
        }
        #endif
        #if 1
        if (ghostInter[])
        {
            double ghostf = 0., totalWeight = 0., alphaf1 = 0;
            coord ghostCell = {x,y,z}, nf2;
            int count = 0;

            // 3.c. extrapolate f from direct neighbors to get initial guess
            int cond0 = fr0[] <= 1e-10;
            //int cond1 = fr0[] >= ibm[] - 1e-10;
            coord ns0 = {-ns.x[], -ns.y[], -ns.z[]};
            foreach_neighbor() {
                if (extra[]) {
                    coord nf1 = {nf.x[], nf.y[], nf.z[]}, ns1 = {-ns.x[], -ns.y[], -ns.z[]};
                    if (cond0 && dot_product_norm(nf1, ns0) < 0)
                        continue;
                    //if (cond1 && dot_product(nf1, ns0) > 0)
                    //    continue;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }
                    
                    double cellWeight = ibm[] * (1. - ibm[]) * f0[] * (1. - f0[]);

                    totalWeight += cellWeight;

                    ghostf += cellWeight * rectangle_fraction (nf1, alphaf[], leftPoint, rightPoint);

                    nf2 = nf1;
                    alphaf1 = alphaf[];
                    #if PRINTCA && 0
                    double product = dot_product_norm(nf1, ns0);
                    double angle = dot_product_angle(nf1, ns0);
                    fprintf(stderr, "extrapolate: (%g, %g, %g) nf={%g, %g, %g} ns={%g, %g, %g} alphaf=%g cw=%g ghostf=%g prod=%g angle=%g cond0=%d cond1=%d\n",
                        x, y, z, nf1.x, nf1.y, nf1.z, ns0.x, ns0.y, ns0.z, alphaf[], cellWeight, ghostf, product, angle, cond0, cond1);
                    #endif

                    count++;
                }
            }

            #if PRINTCA && 0
            fprintf(stderr, " 1st pre (%g, %g, %g) tw = %g gf = %g\n", x, y, z, totalWeight, ghostf);
            #endif

            // 3.d. if extrapolation results in f inside the interface cell, adjust it to preserve fr
            if (totalWeight > 0.) {
                #if PRINTCA && 0
                coord ns1 = {ns.x[], ns.y[], ns.z[]}, nf1 = {nf.x[], nf.y[], nf.z[]};
                fprintf (stderr, "1st %d (%g, %g, %g) nf ={%g,%g} f[]=%0.15g fr0[]=%0.15g tw=%g gf=%g ", 
                                  count, x, y, z, nf2.x, nf2.y, f[], fr0[], totalWeight, ghostf);
                fprintf (stderr, "gf/tw=%g alphaf1=%g\n", ghostf/totalWeight, alphaf1);
                #endif

                f[] = ghostf / totalWeight;
                if (f[] > 1 - INT_TOL) f[] = 1;

                #if PRINTCA
                fprintf (stderr, "  f`[]=%0.15g \n", f[]);
                #endif
            }
            // 3.e otherwise, cell should not be used to enforce C.A.
            else if (ghostf <= 0 && fr0[] <= 1e-10 ) {
                f[] = 0; // was 1 instead of fr0!
            }
        }
        #endif
    }

    boundary({f});
    reconstruction(f,nf,alphaf); // this isn't necessary or even used!

    scalar ftmp[];
    foreach() {
        ftmp[] = f[];
    }
    boundary({ftmp,fr0});

#if PRINTCA
    fprintf(stderr, "Starting full ghost cell correction\n");
#endif
#if 1
    foreach(serial) {
        //gf2[] = f[];
        if (ghostInter[] && fr0[] >= ibm[] - 1e-10 && f[] != f0[] && f[] != 1.) {
            #if 1
            coord mf = interface_normal(point, ftmp);
            if (!mf.x && !mf.y && !mf.z) {
                double ftemp = ftmp[];
                ftmp[] = 0.5;
                mf = interface_normal(point, ftmp);
                ftmp[] = ftemp;
            }
            #endif
            coord ns1 = {ns.x[], ns.y[], ns.z[]};
            double alpha = plane_alpha(ftmp[], mf);           
            double fr = immersed_fraction(ftmp[], mf, alpha, ns1, alphas[], 
                                         (coord){-0.5,-0.5,-0.5}, (coord){0.5,0.5,0.5})*ibm[];
            if (!approx_equal_double(fr, fr0[])) {
                const tripoint tcell = fill_tripoint (fr0[], mf, alpha, ns1, alphas[], f[], ibm[]);
                double alphatmp = ghost_alpha(tcell, -0.6, 0.6);
                f[] = plane_volume(mf, alphatmp);
            }
        #if PRINTCA
            fprintf(stderr, "(%g, %g, %g) f=%g fr=%g ibm=%g n={%g, %g, %g} ns={%g, %g, %g} as=%g\n",
                x, y, z, f[], fr0[], ibm[], mf.x, mf.y, mf.z, ns1.x, ns1.y, ns1.z, alphas[]);
        #endif
        }
    }
    foreach() {
        ftmp[] = f[];
    }
    boundary({f,ftmp});
    #endif
#if PRINTCA
    fprintf(stderr, "Starting three phase cell correction\n");
#endif
#if 1
    foreach(serial) {
        if (extra[]) {
            coord nf1 = interface_normal(point, ftmp);
            if (!nf1.x && !nf1.y && !nf1.z) {
              f[] = 0;
              continue;
            }
            coord ns1 = {ns.x[], ns.y[], ns.z[]};
            double alpha = immersed_alpha_temp (f[], ibm[], nf1, alphaf[], ns1, alphas[], fr0[]);
            f[] = plane_volume(nf1, alpha);
        }
    }
#endif
#if PRINTCA
    fprintf(stderr, "Starting three phase cell correction done\n");
#endif

    boundary ({f});
}

