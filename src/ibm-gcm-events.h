
/*
"vof" event executes at the beginning of the event loop (before advection) but
should execute after the new volume fraction fields have been initalized.

Its main purpose is to update the metric fields to account for a moving interface.
*/

// TODO: this event is very sensitive to MPI and can cause crashes w/AMR
event update_metric (i++)
{
    // update metrics considering immersed boundary

    boundary((scalar *){ibm, ibmf});
    foreach() {
        if (ibm[] > 0.5) // fluid cell
            ibmCells[] = 1;
        else // ghost or solid cell
            ibmCells[] = 0;
    }
    trash ({ibmFaces});
    foreach_face() {
       //if (is_ghost_cell (point, ibm) || ibm[] > 0.5)
        if (ibm[] > 0. || ibm[-1] > 0.)
            ibmFaces.x[] = 1;
        else // solid sell
            ibmFaces.x[] = 0;
    }

    boundary((scalar *){ibmFaces, ibmCells});
}

#if 0
event end_timestep (i++)
{
    // update metrics considering immersed boundary

    boundary((scalar *){ibm, ibmf});
#if 0
    foreach() {
        if (ibm[] > 0.5) // fluid cell
            ibmCells[] = 1;
        else // ghost or solid cell
            ibmCells[] = 0;
    }
    boundary({ibmCells});
#endif
#if 1
    //trash ({ibmFaces});
    foreach_face() {
       //if (is_ghost_cell (point, ibm) || ibm[] > 0.5)
        if (ibm[] > 0. || ibm[-1] > 0.)
            ibmFaces.x[] = 1;
            //ibmFaces.x[] = ibmFaces.x[];
        else // solid sell
            ibmFaces.x[] = 0;
    }

    //fm = ibmFaces;
    //restriction({fm});
    boundary((scalar *){ibmFaces, ibmCells});
#endif
}
#endif

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

vector normals[];
vector midPoints[];

event acceleration (i++)
{

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
   
    // 2. Identify ghost cells and calculate and assign their values to enforce B.C
    foreach() {
        if (is_ghost_cell(point, ibm)) {
           fragment interFrag;
           coord fluidCell, ghostCell = {x,y,z};

           closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell);
           coord boundaryInt = boundary_int (point, interFrag, fluidCell, ibm);
           coord imagePoint = image_point (boundaryInt, ghostCell);
    
           coord imageVelocity = image_velocity (point, u, imagePoint, midPoints);
           double bix = boundaryInt.x, biy = boundaryInt.y, biz = boundaryInt.z;
           foreach_dimension() {
               u.x[] = 2 * uibm_x(bix, biy, biz) - imageVelocity.x;
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

    boundary((scalar *){u, p, pf});
}


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

event end_timestep (i++)
{
    //vector normals[];
    //vector midPoints[];

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

    boundary((scalar *){midPoints, normals});

    correction(-dt);  // remove old pressure from velocity field
    boundary((scalar *){u});

    // 2. Apply the pressure B.C
    foreach() {
        if (is_ghost_cell(point, ibm)) {
           fragment interFrag;
           coord fluidCell, ghostCell = {x,y,z};

           closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell);
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
    
    boundary({p, pf});
    centered_gradient (p, g);
    boundary({g});

    correction(dt);  // add new pressure to velocity field
    boundary((scalar *){u});

    // 3. Update the velocity B.C
    foreach() {
        if (is_ghost_cell(point, ibm)) {
           fragment interFrag;
           coord fluidCell, ghostCell = {x,y,z};

           closest_interface (point, midPoints, ibm, normals, &interFrag, &fluidCell);
           coord boundaryInt = boundary_int (point, interFrag, fluidCell, ibm);
           coord imagePoint = image_point (boundaryInt, ghostCell);
    
           coord imageVelocity = image_velocity (point, u, imagePoint, midPoints);
           double bix = boundaryInt.x, biy = boundaryInt.y, biz = boundaryInt.z;
           foreach_dimension() {
               u.x[] = 2 * uibm_x(bix, biy, biz) - imageVelocity.x;
           }
       }
       else if (ibm[] == 0) {
           foreach_dimension() {
               u.x[] = 0.;
           }
       }
    }
    boundary((scalar *){u});
}

