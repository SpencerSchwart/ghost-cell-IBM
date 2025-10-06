
vector bi[];

/*
"vof" event executes at the beginning of the event loop (before advection) but
should execute after the new volume fraction fields have been initalized.

Its main purpose is to update the metric fields to account for a moving interface.
*/

// TODO: this event is very sensitive to MPI and can cause crashes w/AMR


event update_metric (i++)
{
    // update metrics considering immersed boundary
    boundary({ibm, ibmf});
    foreach() {
        if (ibm[] > 0.5) // fluid cell
            ibmCells[] = 1;
        else // ghost or solid cell
            ibmCells[] = 0;
    }

    foreach_face() {
       if (ibm[] > 0 || ibm[-1] > 0 || is_ghost_cell (point, ibm))
            ibmFaces.x[] = 1;
       else
            ibmFaces.x[] = 0;
    }

#if AXI
    cm_update (cm, ibmCells);
    fm_update (fm, ibmFaces);
    boundary ({fm, cm, ibmFaces, ibmCells});
#else
    boundary({ibmFaces, ibmCells});
#endif
}


/*
###### NOT NECESSARY FOR STATIONARY SOLID ######

The advection_term event is overloaded so the one below is executed first. This
is to handle "fresh" cells created by moving boundary. We do this by interpolating
the fresh cell's velocity (and pressure?) to create a smoother transition from
solid/ghost cell --> fluid cell.

Note the volume fraction field of the previous time step, ibm0, is updated in this event.

TODO: 3D implementation
TODO: better verify if its working well
*/

#ifdef MOVING
scalar fresh[];
event advection_term (i++)
{
    vector normals[];
    vector midPoints[];

    // 1. Initalize fields to hold interface normals and fragment midpoints
    //    TODO: this pass can be improved, if not avoided entirely.
    foreach() {
        coord midPoint, n;
        if (on_interface(ibm)) {
            centroid_point (point, ibm, &midPoint, &n);
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = -n.x;
            }
        }
        else if (ibm[] == 1 && empty_neighbor (point, &midPoint, &n, ibm)) {
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = -n.x;
            }
        }
        else {
            foreach_dimension() {
                midPoints.x[] = 0;
                normals.x[] = 0;
            }
        }
    }

    boundary({midPoints, normals});

    // 2. Find fresh cells and calculate their velocity.
    foreach() {
        if (is_fresh_cell(ibm0, ibm)) {
            fragment interFrag;
            coord freshCell = {x,y}, boundaryInt, n;
            fresh[] = 1;
            ibm_geometry (point, &boundaryInt, &n);
            foreach_dimension()
                boundaryInt.x = freshCell.x + boundaryInt.x * Delta;
            coord imagePoint = fresh_image_point (boundaryInt, freshCell);
            coord imageVelocity = image_velocity (point, u, imagePoint, midPoints);
#if 1
            p[] = image_pressure (point, p, imagePoint);
            pf[] = image_pressure (point, pf, imagePoint);
#endif
           double bix = boundaryInt.x, biy = boundaryInt.y, biz = boundaryInt.z;
            foreach_dimension() {
                u.x[] = (imageVelocity.x + uibm_x(bix,biy,biz)) / 2;
            }
        }
        else if (is_ghost_cell(point, ibm)) {
           fragment interFrag;
           coord fluidCell, ghostCell = {x,y};

           closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell);
           coord boundaryIntercept = boundary_int (point, interFrag, fluidCell, ibm);
           coord imagePoint = image_point (boundaryIntercept, ghostCell);
#if 0  
           if (ibm[] > 0 && ibm0[] <= 0) {
               p[] = image_pressure (point, p, imagePoint);
               pf[] = image_pressure (point, pf, imagePoint);
           }
#endif
        }
        else
            fresh[] = 0;
    }

    foreach() {
        ibm0[] = ibm[];
    }
    boundary({u});
}
#endif


/*
The accleration event is exectued after the advection and viscous events, so u here
is the intermediate velocity field. The purpose of this event is to update the ghost
cell's velocity to better enforce the no slip B.C. since u has changed.

TODO: is calculating and assigning pressure necessary here?
TODO: is assigning pressure to full ghost cells necessary? probably not
*/

scalar ibalphas[];
vector normals[];
vector midPoints[];
//vector bi[];

#if 1
event acceleration (i++)
{
    trash({normals, midPoints});
    // 1. Initalize fields to hold interface normals and fragment midpoints
    //    TODO: this pass can be improved, if not avoided entirely.
    foreach() {
        coord midPoint, n;
        if (on_interface(ibm)) {
            centroid_point (point, ibm, &midPoint, &n, &ibalphas[]);
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = -n.x;
            }
            ibalphas[] *= -1;
        }
        else if (ibm[] == 1 && empty_neighbor (point, &midPoint, &n, ibm)) {
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = -n.x;
            }
            ibalphas[] = 0;
        }
        else {
            foreach_dimension() {
                midPoints.x[] = 0;
                normals.x[] = 0;
            }
            ibalphas[] = 0;
        }
    }

    boundary({midPoints, normals, ibalphas});
  
     // ghost boundary intercepts
    u.x.mp = bi; u.y.mp = bi; 
#if dimension == 3
    u.z.mp = bi;
#endif

    // 2. Identify ghost cells and calculate and assign their values to enforce B.C
    foreach() {
        if (is_ghost_cell(point, ibm)) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y,z};
            PointIBM bioff = {0,0,0};

            closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell, &bioff);
            coord boundaryInt = boundary_int (point, interFrag, fluidCell, ibm);

            coord imagePoint = image_point (boundaryInt, ghostCell);

            foreach_dimension()
                bi.x[] = boundaryInt.x;

            coord imageVelocity = image_velocity2 (point, u, imagePoint, bioff, 
                                                   midPoints, normals, ibalphas);
            
            if (local_bc_coordinates) {
                coord n = interFrag.n, t1, t2;
                normal_and_tangents (&n, &t1, &t2);

                coord projVelocity = {dot_product(imageVelocity, n),
                                      dot_product(imageVelocity, t1),
                                      dot_product(imageVelocity, t2)};

                coord gcProjVelocity = {0,0,0};
                foreach_dimension() {
                    bool bctype[2] = {false, false};
                    double vb = u.x.boundary[immersed] (point, point, u.x, bctype);
                    bool dirichlet = bctype[0], nslip = bctype[1];
                    if (dirichlet) {
                        if (nslip) {
                            double delta = 0;
                            foreach_dimension()
                                delta += sq(ghostCell.x - imagePoint.x);
                            delta = sqrt(delta);

                            // only for stationary solids right now (hence 0 -)
                            gcProjVelocity.x = delta*(0 - projVelocity.x)/(0.5*delta + vb) + 
                                               projVelocity.x;
                        }
                        else
                            gcProjVelocity.x = 2*vb - projVelocity.x;
                    }
                    else {
                        gcProjVelocity.x = projVelocity.x;
                    }
                }

                double gcn = gcProjVelocity.x, gct1 = gcProjVelocity.y, gct2 = gcProjVelocity.z;
                foreach_dimension()
                    u.x[] = gcn*n.x + gct1*t1.x + gct2*t2.x;
            }
            else { // !local_bc_coordinates, i.e. use n ≡ x, t ≡ y, and r ≡ z.
                foreach_dimension() {
                    bool dirichlet = false;
                    double vb = u.x.boundary[immersed] (point, point, u.x, &dirichlet);
                    if (dirichlet) {
                        //fprintf(stderr, "WARNING: Navier boundary condtion only available for local bc coordinates\n");
                        u.x[] = 2*vb - imageVelocity.x;
                    }
                    else {
                        u.x[] = imageVelocity.x;
                    }
                }
            }
            if (ibm[] <= 0.) { // is pressure b.c. necessary here?
                p[] = image_pressure (point, p, imagePoint);
                pf[] = image_pressure (point, pf, imagePoint);
            }
       }
       else if (ibm[] == 0) {
           p[] = 0; pf[] = 0;
           foreach_dimension() {
               u.x[] = 0.;
           }
       }
    }

    boundary((scalar *){u, p, pf, bi});
}
#endif

/*
end_timestep updates the boundary conditions (both pressure and velocity) since
both fields are altered after projection and velocity correction.

Note: we only explicity apply pressure B.C to ghost cells which are entirely filled
with solid (ibm = 0) since the pressure solver can't set them. This means that ALL
cells with a fragment of interface are assigned their pressure via the projection step.

TODO: make it so midPoints field doesn't need to be recalculated. Replace it
      with ad-hoc offset function.
TODO: are all of these boundary()'s necessary?
TODO: is assigning pressure to full ghost cells necessary?
*/

#if 1
event end_timestep (i++)
{
    correction(-dt);  // remove old pressure from velocity field

    // 2. Apply the pressure B.C
    foreach() {
        if (is_ghost_cell(point, ibm)) {
           fragment interFrag;
           coord fluidCell, ghostCell = {x,y,z};
           PointIBM bioff;

           closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell, &bioff);
           coord boundaryInt = boundary_int (point, interFrag, fluidCell, ibm);
           coord imagePoint = image_point (boundaryInt, ghostCell);
    
           if (ibm[] <= 0.) {
               p[] = image_pressure (point, p, imagePoint);
               pf[] = image_pressure (point, pf, imagePoint);
           }
       }
       else if (ibm[] == 0) {
           p[] = 0; pf[] = 0;
       }
    }
    
    boundary({u, p, pf});
    centered_gradient (p, g);
    //boundary({g}); // is this necessary?

    correction(dt);  // add new pressure to velocity field
    boundary({u});

    u.x.mp = bi; u.y.mp = bi; 
#if dimension == 3
    u.z.mp = bi;
#endif

    // 3. Update the velocity B.C
    foreach() {
        if (is_ghost_cell(point, ibm)) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y,z};
            PointIBM bioff;

            closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell, &bioff);
            coord boundaryInt = boundary_int (point, interFrag, fluidCell, ibm);
            coord imagePoint = image_point (boundaryInt, ghostCell);
  
            foreach_dimension()
                bi.x[] = boundaryInt.x;

            coord imageVelocity = image_velocity2 (point, u, imagePoint, bioff, 
                                                   midPoints, normals, ibalphas);

            if (local_bc_coordinates) {
                coord n = interFrag.n, t1, t2;
                normal_and_tangents (&n, &t1, &t2);

                coord projVelocity = {dot_product(imageVelocity, n),
                                      dot_product(imageVelocity, t1),
                                      dot_product(imageVelocity, t2)};
                coord gcProjVelocity = {0,0,0};
                foreach_dimension() {
                    bool bctype[2] = {false, false};
                    double vb = u.x.boundary[immersed] (point, point, u.x, bctype);
                    bool dirichlet = bctype[0], nslip = bctype[1];
                    if (dirichlet) {
                        if (nslip) {
                            double delta = 0;
                            foreach_dimension()
                                delta += sq(ghostCell.x - imagePoint.x);
                            //delta = sqrt(delta)/Delta; // normalize by cell size
                            delta = sqrt(delta);

                            // only for stationary solids right now (hence 0 -)
                            gcProjVelocity.x = delta*(0 - projVelocity.x)/(0.5*delta + vb) + 
                                               projVelocity.x;
                        }
                        else
                            gcProjVelocity.x = 2*vb - projVelocity.x;
                    }
                    else { // neumann
                        gcProjVelocity.x = projVelocity.x;
                    }
                }

                double gcn = gcProjVelocity.x, gct1 = gcProjVelocity.y, gct2 = gcProjVelocity.z;
                foreach_dimension()
                    u.x[] = gcn*n.x + gct1*t1.x + gct2*t2.x;
            }
            else { // !local_bc_coordinates, i.e. use n ≡ x, t ≡ y, and r ≡ z.
                foreach_dimension() {
                    bool dirichlet = false;
                    double vb = u.x.boundary[immersed] (point, point, u.x, &dirichlet);
                    if (dirichlet) {
                        u.x[] = 2*vb - imageVelocity.x;
                    }
                    else {
                        u.x[] = imageVelocity.x;
                    }
                }
            }
       }
       else if (ibm[] == 0) {
           foreach_dimension() {
               u.x[] = 0.;
           }
       }
    }
    boundary((scalar *){u, bi});
}
#endif
