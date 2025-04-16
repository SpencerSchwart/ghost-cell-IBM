
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
    if (dimension == 2) {
        if (- ns.x * nf.y + ns.y * nf.x > 0) { // 2D cross product
            n.x = - ns.x * cos(angle) + ns.y * sin(angle);
            n.y = - ns.x * sin(angle) - ns.y * cos(angle);
        }
        else {
            n.x = - ns.x * cos(angle) - ns.y * sin(angle);
            n.y =   ns.x * sin(angle) - ns.y * cos(angle);
        }
    }
    else { // dimension == 3
    #if 1 // spherical coordinates
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
        coord rot_axis;     //Quaternion and rotation matrix 

        foreach_dimension() //rotation axis to define a quaternion
           rot_axis.x = ns.y*nf.z - ns.z*nf.y; 
         normalize(&rot_axis);

         coord p1 = {ns.x, ns.y, ns.z};  //original position of point
    
         double re = cos(angle/2.);           //real part of quaternion
         double x = rot_axis.x*sin(angle/2.); //imaginary i part of quaternion
         double y = rot_axis.y*sin(angle/2.); //imaginary j part of quaternion
         double z = rot_axis.z*sin(angle/2.); //imaginary k part of quaternion

         n.x = re*re*p1.x + 2*y*re*p1.z - 2*z*re*p1.y + x*x*p1.x + 2*y*x*p1.y + 2*z*x*p1.z - z*z*p1.x - y*y*p1.x;
         n.y = 2*x*y*p1.x + y*y*p1.y + 2*z*y*p1.z + 2*re*z*p1.x - z*z*p1.y + re*re*p1.y - 2*x*re*p1.z - x*x*p1.y;
         n.z = 2*x*z*p1.x + 2*y*z*p1.y + z*z*p1.z - 2*re*y*p1.x - y*y*p1.z + 2*re*x*p1.y - x*x*p1.z + re*re*p1.z;
    #endif
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
            double mag = fabs(nc.x) + fabs(nc.y) + fabs(nc.z) + SEPS;
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
            if (!fluidNeighbors) 
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
    vector n[];
    scalar alpha[];

    // 1. reconstruct n and alpha considering the contact angle B.C.
    reconstruction_contact (f, n, alpha);

    // 2. look for pure solid cells near triple point to enforce B.C
    foreach() {
        if (ibm[] <= 0.) {  // pure solid cell
            double ghostf = 0., totalWeight = 0.; // ghostf = volume fraction of ghost cell
            coord ghostCell = {x,y,z};
    
            // 2.a. Calculate weights
            foreach_neighbor() {
                if (on_interface(ibm) && on_interface(f)) {

                    double cellWeight = ibm[] * (1. - ibm[]) * f[] * (1. - f[]);
                    totalWeight += cellWeight;

                    coord leftPoint = {x, y, z}, rightPoint;
                    foreach_dimension() {
                        leftPoint.x = (ghostCell.x - leftPoint.x) / Delta - 0.5;
                        rightPoint.x = leftPoint.x + 1.;
                    }

                    coord nf = {n.x[], n.y[], n.z[]};
                    ghostf += cellWeight * rectangle_fraction (nf, alpha[], leftPoint, rightPoint);
                }
            }
            
            // 2.b. Assign volume fraction in ghost/solid cell
            if (totalWeight > 0.)
                f[] = ghostf / totalWeight;
        }
    }

    //3. remove unnecessary fluid
    clean_fluid(f, ibm);
    boundary ({f});
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
