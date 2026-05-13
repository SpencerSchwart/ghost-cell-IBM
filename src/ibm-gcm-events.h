
vector bis[];

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
       //if (cs[] > 0 || cs[-1] > 0)
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

scalar ibalphas[];
vector normals[];
vector midPoints[];

coord interpolate_image_point(coord ip, coord nsg, int * rank = NULL) 
{
    //double ux = nodata, uy=nodata;
    coord uip = {nodata};
    //foreach_point(ip.x, ip.y, reduction(min:ux) reduction(min:uy)) {
    foreach_point_ibm(ip.x, ip.y, ip.z, rank, serial) {

        //fprintf(fp, "(%g, %g): point = (%g, %g) pid=%g %g\n", x, y, ip0.x, ip0.y, pid(), pid[]);
        //fflush(fp);
        //foreach_dimension()
        //    bi.x[] = bis[i].x;

        uip = image_velocity (point, u, ip, nsg, midPoints, normals, ibalphas);

        //uip[i] = image_velocity_noff (point, u, ips[i], bff[i], nsg[i], midPoints, normals, ibalphas);
	//ux = uip.x, uy = uip.y;
    }
    //return (coord){ux,uy};
    return uip;
}

void fill_interface_data() {
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
}

scalar gid[], gbid[];
scalar pid[];

#if _MPI

typedef struct {
    int id;    // ID in array
    double x, y, z;
    double nx, ny, nz;
} MPIcoord;

#endif

void update_gc_velocity()
{

    foreach()
      pid[] = pid();

    boundary({pid});

    u.x.mp = bis; u.y.mp = bis; 
#if dimension == 3
    u.z.mp = bis;
#endif

#if 1

    /**
    First count all ghost cells. */
    int gcount = 0, gbcount = 0;

    Array * gcid  = array_new();

    foreach(serial, reduction(+:gcount) reduction(+:gbcount)) {
        gid[] = gbid[] = -1;
        if (is_ghost_cell(point, cs)) {
            gid[] = gcount;
            array_append(gcid, &gcount, sizeof(int));
            gcount++;
        }
    }

#if 0
    char name[80];
    sprintf(name, "%d-out", pid());
    FILE * fp = fopen(name, "w");
#endif

    /**
    Allocate arrays for storing the, IP, IP velocity, and normals for each GC */
    long nl  = gcid->len/sizeof(int);
#if 1

    coord ips[nl];
    coord uip[nl];
    coord nsg[nl];

    foreach() {
        if (gid[] > -1) {
            fragment interFrag;
            coord fc, gc = {x,y,z};
            PointIBM bioff;

            coord bi, ip, n;
            if (cs.ibm.bi) 
                bi = cs.ibm.bi(gc);
            else {
                closest_interface (point, midPoints, cs, normals, ibalphas, &interFrag, &fc, &bioff);
                bi = boundary_int (point, interFrag, fc, cs);
                n = interpolate_normal (point, bi, fc, bioff, midPoints, normals);
            }
            if (cs.ibm.normal)
                n = cs.ibm.normal(gc);
            ip = image_point (bi, gc, n);

            nsg[(int)gid[]] = n;
            ips[(int)gid[]] = ip;
        }
    }

#endif

    /**
    When using MPI, an image point may lay in another rank, requiring extra communication. 
    We identify these problematic points and store its index in the bips (bad image points) array.
    We also have to store the rank the IP belongs to in an array, bpid. */

#if _MPI
    Array * bips  = array_new();
    Array * bpid  = array_new();
#endif

#if 1
    for (int i = 0; i < nl; i++) {
	//fprintf(fp, "(%g, %g) pid=%g", ips[i].x, ips[i].y, pid());
	//fflush(fp);

        int rank = -1;
        uip[i] = interpolate_image_point(ips[i], nsg[i], &rank);

        if (uip[i].x == nodata) {
            fprintf(fp, "HELP! GC rank = %d, IP rank = %d\n", pid(), rank);
            array_append(bips, &i, sizeof(int));
            array_append(bpid, &rank, sizeof(int));
        }

        fflush(fp);

	//fprintf(fp, " uip=(%g, %g)\n", uip[i].x, uip[i].y);
#if 0
        coord ip0 = ips[i];
        foreach_point(ip0.x, ip0.y, reduction(min:val)) {
	    fprintf(fp, "(%g, %g): point = (%g, %g) pid=%g %g\n", x, y, ip0.x, ip0.y, pid(), pid[]);
	    fflush(fp);
            foreach_dimension()
                bis.x[] = bi[i].x;
            //uip[i] = image_velocity (point, u, ips[i], nsg[i], midPoints, normals, ibalphas);
            //uip[i] = image_velocity_noff (point, u, ips[i], bff[i], nsg[i], midPoints, normals, ibalphas);
        }
#endif
    }

    /**
    If "bad" image points exists, we count the total number for each rank and store it
    in an array whose indices correspond to all active processes. We then communicate this
    array to processors. */

#if _MPI
    int scount[npe()]; // send
    int rcount[npe()]; // receive

    for (int r = 0; r < npe(); r++) {
        scount[r] = 0;
        rcount[r] = 0;
    }

    int send_total = bips->len/sizeof(int);

    int * ip_rank = bpid->p;
    int * ip_ids = bips->p;

    for (int i = 0; i < send_total; i++) {
        scount[ip_rank[i]]++;
    }

    MPI_Alltoall(scount, 1, MPI_INT, rcount, 1, MPI_INT, MPI_COMM_WORLD);

    /**
    Calcuate the total number of receiving indices. Note we already calculated
    the number of ones being sent. Then we allocate an array for points we will send
    and receive. */

    int recv_total = 0;
    for (int r = 0; r < npe(); r++) 
        recv_total += rcount[r];

    MPIcoord * send_ips = malloc (send_total * sizeof(MPIcoord));
    MPIcoord * recv_ips = malloc (recv_total * sizeof(MPIcoord));

    /**
    After we exchange how many "bad" IPs each process both has (send) and needs to interpolate
    for other processes (receive), we need to calculate the starting index for each process in
    the count arrays. We store these starting index in these displacement arrays. */

    int sdispls[npe()];
    int rdispls[npe()];

    /**
    The starting index for pid 0 is always 0. */

    sdispls[0] = 0;
    rdispls[0] = 0;

    for (int r = 1; r < npe(); r++) {
        sdispls[r] = sdispls[r - 1] + scount[r - 1];
        rdispls[r] = rdispls[r - 1] + rcount[r - 1];
    }

    /**
    We now need to fill the send_ips array with all of the image point coordinates.
    However, we need to make sure it is sorted and grouped by rank so that it
    corresponds to the displacement array sdispls. To do so, we create a temporary 
    copy. */

    int sdisp_temp[npe()];
    
    for (int r = 0; r < npe(); r++)
        sdisp_temp[r] = sdispls[r];

    for (int i = 0; i < send_total; i++) {
        int id = ip_ids[i];  // grab id in main array
        coord ipt = ips[id]; // grab IP
        coord n = nsg[id];   // grab normal

        int rank = ip_rank[i];  // grab IP home rank
        
        int index = sdisp_temp[rank];
        send_ips[index] = (MPIcoord){id, ipt.x, ipt.y, ipt.z, n.x, n.y, n.z};

        sdisp_temp[rank]++;
    }

    /**
    Before we can communicate the IP to their respective rank, we must first
    convert the count and displacement arrays into byte format since user defined
    types, like MPIcoord, are not natively supported in MPI. */

    int scount_byte[npe()];
    int rcount_byte[npe()]; 
    int sdispls_byte[npe()];
    int rdispls_byte[npe()];

    for (int r = 0; r < npe(); r++) {
        scount_byte[r] = scount[r] * sizeof(MPIcoord); 
        rcount_byte[r] = rcount[r] * sizeof(MPIcoord); 
        sdispls_byte[r] = sdispls[r] * sizeof(MPIcoord);
        rdispls_byte[r] = rdispls[r] * sizeof(MPIcoord);
    }

    MPI_Alltoallv(send_ips, scount_byte, sdispls_byte, MPI_BYTE,
                  recv_ips, rcount_byte, rdispls_byte, MPI_BYTE,
                  MPI_COMM_WORLD);

    /**
    Interpolate all of the received IPs and store them in an array to be communicated.
    Note we do not have to sort this array, since recv_ips is already sorted. */

    MPIcoord * uip_send = malloc(recv_total * sizeof(MPIcoord));

    for (int i = 0; i < recv_total; i++) {

        MPIcoord data = recv_ips[i];
        coord ip = {data.x, data.y, data.z};
        coord n = {data.nx, data.ny, data.nz};
        
        int rank = -1;
        coord uip = interpolate_image_point(ip, n, &rank);

        assert(rank == pid());

        uip_send[i] = (MPIcoord){data.id, uip.x, uip.y, uip.z};
    }

    MPIcoord * uip_recv = malloc(send_total * sizeof(MPIcoord));

    MPI_Alltoallv(uip_send, rcount_byte, rdispls_byte, MPI_BYTE,
                  uip_recv, scount_byte, sdispls_byte, MPI_BYTE,
                  MPI_COMM_WORLD);

    /**
    The interpolated IP values originally requested are now stored in uip_recv.
    We store them back in the main uip array with the rest of the data. */

    for (int i = 0; i < send_total; i++) {
        MPIcoord data = uip_recv[i];
        uip[data.id] = (coord){data.x, data.y, data.z};
    }

#endif

    fflush(fp);
    fclose(fp);
#endif

    array_free(gcid);

#if _MPI
    if (uip_send) free(uip_send);
    if (uip_recv) free(uip_recv);
    if (send_ips) free(send_ips);
    if (recv_ips) free(recv_ips);
    
    array_free(bips);
    array_free(bpid);
#endif

#endif

    scalar gcin[];

    //bool avg_deep = !cs.ibm.bi && !cs.ibm.expression && !cs.ibm.normal;
    bool avg_deep = false;
    //bool avg_deep = true;

    foreach() {
        gcin[] = -1;
        if (is_ghost_cell(point, cs)) {

            bool inside = true;
            foreach_direct_neighbor() 
                if (gc[])
                    inside = false;

            gcin[] = inside? 1: 0;
            
            if (avg_deep && inside) {
                foreach_dimension()
                    u.x[] = 0;
                continue;
            }

            fragment interFrag;
            coord fc, gc = {x,y,z};
            PointIBM bioff;

            coord bi = {0,0,0}, ip, n;
            if (cs.ibm.bi)
                bi = cs.ibm.bi(gc);
            else {
                closest_interface (point, midPoints, cs, normals, ibalphas, &interFrag, &fc, &bioff);
                bi = boundary_int (point, interFrag, fc, cs);

                n = interpolate_normal (point, bi, fc, bioff, midPoints, normals);
            }
            if (cs.ibm.normal)
                n = cs.ibm.normal(gc);
            ip = image_point(bi, gc, n);

            foreach_dimension()
                bis.x[] = bi.x;

            //coord imageVelocity = image_velocity (point, u, ip, n, 
            //                                      midPoints, normals, ibalphas);

            coord imageVelocity = uip[(int)gid[]];

            if (local_bc_coordinates) {
                coord t1, t2;
                normal_and_tangents (&n, &t1, &t2);

                coord projVelocity = {dot_product(imageVelocity, n),
                                      dot_product(imageVelocity, t1),
                                      dot_product(imageVelocity, t2)};

                coord gcProjVelocity = {0,0,0};

                foreach_dimension() {
                    bool bctype[2] = {false, false};
                    double vb = u.x.boundary[immersed] (point, point, u.x, bctype);
                    bool dirichlet = bctype[0], navier = bctype[1];
                    if (dirichlet) {
                        if (navier) {
                            double delta = 0;
                            foreach_dimension()
                                delta += sq(gc.x - ip.x);
                            delta = sqrt(delta)/(2*Delta);
                            vb /= Delta;

                            gcProjVelocity.x = -projVelocity.x*(delta - vb)/(delta + vb);
                        }
                        else {
                            gcProjVelocity.x = 2*vb - projVelocity.x;
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

#if 0
    if (avg_deep) {
        foreach() {
            if (gcin[] == 1) {

                coord ugc = {0,0};
                int count = 0;
                foreach_near_neighbor() {
                   if (u.x[]) {
                       foreach_dimension()
                           ugc.x += u.x[];
                       count++;
                   }
               }
               if (count) {
                   foreach_dimension()
                       u.x[] = ugc.x/(double)count;
               }
            }
        }
    }
#endif
    boundary((scalar *){u, bis});
}

event end_timestep(i++)
{
    fill_interface_data(); 
    update_gc_velocity();
}

#if 0
event acceleration (i++)
{
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

#if 0
    // 2. Identify ghost cells and calculate and assign their values to enforce B.C
    foreach() {
        if (is_ghost_cell(point, cs)) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y,z};
            PointIBM bioff = {0,0,0};

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

#if 0

scalar gid[];

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

#if 0
    // Identify and mark ghost cells
    int gcount = 0;

    Array * gcid = array_new();

    foreach(serial, reduction(+:gcount)) {
        gid[] = -1;
        if (is_ghost_cell(point, cs)) {
            gid[] = gcount;
            array_append(gcid, &gcount, sizeof(int));
            gcount++;
        }
    }

    char name[80];
    sprintf(name, "%d-out", pid());
    FILE * fp = fopen(name, "w");

    long nl = gcid->len/sizeof(int);
    long offset = 0;

#if _MPI
    MPI_Exscan(&nl, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (pid() == 0)
        offset = 0;
#endif

    int * ap = gcid->p;

#if 0
    foreach(serial) {
        if (gid[] > -1) {
            gid[] += offset;
        }
    }

    for (long i = 0; i < nl; i++) {
        int old = ap[i];
        ap[i] += offset;
        fprintf(fp, "%d nl=%d gcid_old = %d gcid[%d] = %d offset = %d\n", pid(), nl, old, i, ap[i], offset);
    }
    fclose(fp);
#endif

    coord ips[nl];
    coord uip[nl];
    PointIBM bff[nl];
    coord nsg[nl];
    coord bis[nl];

    foreach() {
        if (gid[] > -1) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y,z};
            PointIBM bioff;

            closest_interface (point, midPoints, cs, normals, &interFrag, &fluidCell, &bioff);
            coord boundaryInt = boundary_int (point, interFrag, fluidCell, cs);

            bis[(int)gid[]] = (coord)boundaryInt;
            nsg[(int)gid[]] = (coord)interFrag.n;
            bff[(int)gid[]] = (PointIBM)bioff;
            ips[(int)gid[]] = image_point (boundaryInt, ghostCell);
        }
    }

#if 0
    for (long i = 0; i < nl; i++) {
        coord ip0 = ips[i];
        foreach_point(serial, ip0.x, ip0.y) {
            foreach_dimension()
                bi.x[] = bis[i].x;
            uip[i] = image_velocity_noff (point, u, ips[i], bff[i], nsg[i], midPoints, normals, ibalphas);
            fprintf(fp, "%d (%g, %g) nl=%d gcid[%d]=%d ips={%g, %g} uip={%g, %g} offset=%d %g\n", pid(), x, y, nl, i, ap[i], ips[i].x, ips[i].y, uip[i].x, uip[i].y, offset, cs[]);
        }
    }
    fflush(fp);
    fclose(fp);
#endif

    array_free(gcid);

#if 0
    // should move them to heap allocation using malloc
    //coord bip[gcount]; // boundary intercept points
    //coord uip[gcount]; // image point velocities
    
    //
    foreach() {
        if (gid[] > -1) {
            MPI_Allreduce (MPI_IN_PLACE, size, n, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

            MPI_Exscan (
            gid[] = 
        }
    }
#endif
    
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
#if 0                                                  
            coord imageVelocity0 = uip[(int)gid[]];
            coord imagePoint0 = ips[(int)gid[]];
            PointIBM bioff0 = bff[(int)gid[]];
            
            if (!approx_equal_double(imageVelocity.x,imageVelocity0.x) || !approx_equal_double(imageVelocity.y, imageVelocity0.y))
                fprintf(stderr, "%g (%g, %g) uip0 = {%g, %g} uip1 = {%g, %g} cs=%g ip0={%g, %g} ip1={%g, %g} bioff0={%d, %d} bioff1={%d, %d}\n", 
                    gid[], x, y, imageVelocity0.x, imageVelocity0.y, imageVelocity.x, imageVelocity.y, cs[], 
                    imagePoint0.x, imagePoint0.y, imagePoint.x, imagePoint.y, bioff0.i, bioff0.j, bioff.i, bioff.j);

            imageVelocity = imageVelocity0;
#endif

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
                            gcProjVelocity.x = 2*vb - projVelocity.x;
                        }
                    }
                    else { // neumann
                        gcProjVelocity.x = projVelocity.x;
                    }
                }

#if 0
                //double r = sqrt(sq(imagePoint.x) + sq(imagePoint.y));
                double r = sqrt(sq(ghostCell.x) + sq(ghostCell.y));


                if (r > 0.5*(0.5+0.25)) 
                {
                  gcProjVelocity.x = 0;
                  gcProjVelocity.y = -(r*(sq(0.5/r) - 1.)/(sq(0.5/0.25) - 1.));
                  //gcProjVelocity.y = -sq(0.25)/(sq(0.25) + sq(0.5))*(sq(0.5)/r + r);
                }
                else {
                  gcProjVelocity.x = 0;
                  gcProjVelocity.y = (r*(sq(0.5/r) - 1.)/(sq(0.5/0.25) - 1.));
                  //gcProjVelocity.y = sq(0.25)/(sq(0.25) + sq(0.5))*(sq(0.5)/r + r);
                }
#endif
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

