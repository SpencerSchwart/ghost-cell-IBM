
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
    //boundary({cs, fs});
    cs.dirty = true;
    foreach() {
        if (cs[] > GCV) // fluid cell
            gc[] = 1;
        else // ghost or solid cell
            gc[] = 0;
    }

    foreach_face() {
       if (cs[] > 0 || cs[-1] > 0 || is_ghost_cell (point, cs))
            gcf.x[] = 1;
       else
            gcf.x[] = 0;
    }

#if AXI
    cm_update (cm, gc);
    fm_update (fm, gcf);
    boundary ({fm, cm, gcf, gc});
#else
    boundary({gcf, gc});
#endif
}


#if 0
/*
###### NOT NECESSARY FOR STATIONARY SOLID ######

The advection_term event is overloaded so the one below is executed first. This
is to handle "fresh" cells created by moving boundary. We do this by interpolating
the fresh cell's velocity (and pressure?) to create a smoother transition from
solid/ghost cell --> fluid cell.

Note the volume fraction field of the previous time step, cs0, is updated in this event.

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
        if (on_interface(cs)) {
            centroid_point (point, cs, &midPoint, &n);
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = -n.x;
            }
        }
        else if (cs[] == 1 && empty_neighbor (point, &midPoint, &n, cs)) {
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
        if (is_fresh_cell(cs0, cs)) {
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
        else if (is_ghost_cell(point, cs)) {
           fragment interFrag;
           coord fluidCell, ghostCell = {x,y};

           closest_interface (point, midPoints, cs, normals, &interFrag, &fluidCell);
           coord boundaryIntercept = boundary_int (point, interFrag, fluidCell, cs);
           coord imagePoint = image_point (boundaryIntercept, ghostCell);
#if 0  
           if (cs[] > 0 && cs0[] <= 0) {
               p[] = image_pressure (point, p, imagePoint);
               pf[] = image_pressure (point, pf, imagePoint);
           }
#endif
        }
        else
            fresh[] = 0;
    }

    foreach() {
        cs0[] = cs[];
    }
    boundary({u});
}
#endif
#endif

/*
The accleration event is exectued after the advection and viscous events, so u here
is the intermediate velocity field. The purpose of this event is to update the ghost
cell's velocity to better enforce the no slip B.C. since u has changed.

TODO: is calculating and assigning pressure necessary here?
TODO: is assigning pressure to full ghost cells necessary? probably not
*/

const coord ind = {0,1,2};

scalar ibalphas[];
vector normals[];
vector midPoints[];
#if 1
event acceleration (i++)
//event projection (i++)
{
    trash({normals, midPoints});
    // 1. Initalize fields to hold interface normals and fragment midpoints
    //    TODO: this pass can be improved, if not avoided entirely.
    foreach() {
        coord midPoint, n;
        if (on_interface(cs)) {
            centroid_point (point, cs, &midPoint, &n, &ibalphas[]);
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = -n.x;
            }
            ibalphas[] *= -1;
        }
        else if (cs[] == 1 && empty_neighbor (point, &midPoint, &n, cs)) {
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

#if 0
    vector u2[];
    scalar_clone(u2.x, u.x);
    scalar_clone(u2.y, u.y);

    vector g0[];
    centered_gradient (p, g0);

    foreach() {
        foreach_dimension()
            u2.x[] = u.x[] + gc[]*dt*g0.x[];
    }
#endif

#if 1
    // 2. Identify ghost cells and calculate and assign their values to enforce B.C
    foreach() {
        if (is_ghost_cell(point, cs)) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y,z};
            PointIBM bioff = {0,0,0};

            closest_interface (point, midPoints, cs, normals, &interFrag, &fluidCell, &bioff);
            coord boundaryInt = boundary_int (point, interFrag, fluidCell, cs);
            //rho[] = rho[0,1];

            coord imagePoint = image_point (boundaryInt, ghostCell);
  
            foreach_dimension()
                bi.x[] = boundaryInt.x;

            coord imageVelocity = image_velocity (point, u, imagePoint, bioff, 
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
                            delta = sqrt(delta)/(2*Delta);
                            vb /= Delta;

                            // only for stationary solids right now (hence 0 -)
                            //gcProjVelocity.x = delta*(0 - projVelocity.x)/(0.5*delta + vb) + 
                            //                   projVelocity.x;
                            gcProjVelocity.x = -projVelocity.x*(delta - vb)/(delta + vb);
                        }
                        else {
                            double delta = 0;
                            foreach_dimension()
                                delta += sq(ghostCell.x - imagePoint.x);
                            delta = sqrt(delta)/(2*Delta);

                            double multiplier = delta <= 1? 1: 1;
                            
                            gcProjVelocity.x = 2*vb - projVelocity.x*multiplier;
                        }
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
                        u.x[] = 2*vb - imageVelocity.x;
                    }
                    else {
                        u.x[] = imageVelocity.x;
                    }
                }
            }
            #if 0
            if (cs[] <= 0.) { // is pressure b.c. necessary here?
                p[] = image_pressure (point, p, imagePoint);
                pf[] = image_pressure (point, pf, imagePoint);
            }
            #endif
       }
       else if (cs[] == 0) {
           p[] = 0; pf[] = 0;
           foreach_dimension() {
               u.x[] = 0.;
           }
       }
    }
#endif
    boundary((scalar *){u, p, pf, bi});
}
#endif

/*
end_timestep updates the boundary conditions (both pressure and velocity) since
both fields are altered after projection and velocity correction.

Note: we only explicity apply pressure B.C to ghost cells which are entirely filled
with solid (cs = 0) since the pressure solver can't set them. This means that ALL
cells with a fragment of interface are assigned their pressure via the projection step.

TODO: make it so midPoints field doesn't need to be recalculated. Replace it
      with ad-hoc offset function.
TODO: are all of these boundary()'s necessary?
TODO: is assigning pressure to full ghost cells necessary?
*/

#if 1
event end_timestep (i++)
{
    u.x.mp = bi; u.y.mp = bi; 
#if dimension == 3
    u.z.mp = bi;
#endif

#if 0
    foreach_face() {
        double metric = !fs.x[]? 0: fm.x[]*alpha.x[]/fs.x[];
        uf.x[] -= dt*metric*face_gradient_x (p, 0);
    }
#endif

    // 3. Update the velocity B.C
    foreach() {
        if (is_ghost_cell(point, cs)) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y,z};
            PointIBM bioff;

            closest_interface (point, midPoints, cs, normals, &interFrag, &fluidCell, &bioff);
            coord boundaryInt = boundary_int (point, interFrag, fluidCell, cs);
            coord imagePoint = image_point (boundaryInt, ghostCell);
  
            foreach_dimension()
                bi.x[] = boundaryInt.x;

            coord imageVelocity = image_velocity (point, u, imagePoint, bioff, 
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
                            delta = sqrt(delta)/(2*Delta);
                            vb /= Delta;

                            gcProjVelocity.x = -projVelocity.x*(delta - vb)/(delta + vb);
                        }
                        else {
                             double delta = 0;
                            foreach_dimension()
                                delta += sq(ghostCell.x - imagePoint.x);
                            delta = sqrt(delta)/(2*Delta);

                            double multiplier = delta <= 1? 1: 1;
                            //if (ind.x == 0) multiplier = 1; 
                            //fprintf(stderr, "multiplier = %g delta = %g\n", multiplier, delta);
                            gcProjVelocity.x = 2*vb - projVelocity.x*multiplier;
                        }
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
       else if (cs[] == 0) {
           foreach_dimension() {
               u.x[] = 0.;
           }
       }
    }
    boundary((scalar *){u, bi});
}
#endif
