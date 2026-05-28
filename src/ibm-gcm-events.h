
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

scalar ibalphas[];
vector normals[];
vector midPoints[];

void fill_interface_data() {
    foreach() {
        coord midPoint, n;
        if (cs[] > 0 && cs[] < 1) {
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

scalar gid[], gidd[], gidd0[], dg0[];
scalar pid[];

#if _MPI
typedef struct {
    int id;    // ID in array
    double x, y, z;
    double nx, ny, nz;
} MPIcoord;

#endif

// TODO: doesn't work in openmp?
coord interpolate_image_point (coord ip, coord nsg, int * rank = NULL) 
{
    coord uip = {nodata};
    foreach_image_point (ip.x, ip.y, ip.z, rank)
        uip = image_velocity (point, u, ip, nsg, midPoints, normals, ibalphas);
    return uip;
}

scalar neg[];

void update_gc_velocity()
{
#if 0
    foreach()
      pid[] = pid();
#endif

    /**
    First count all ghost cells. */
    int gcount = 0;

    Array * gcid  = array_new();

    foreach(serial, reduction(+:gcount)) {
        gid[] = gidd0[] = -1;
        gidd[] = 0;
        if (is_ghost_cell(point, cs)) {
            gid[] = gcount;
            gidd[] = (int)is_deep_ghost_cell(point);
            array_append(gcid, &gcount, sizeof(int));
            gcount++;
        }
    }

    /**
    Allocate arrays for storing the, IP, IP velocity, and normals for each GC */

    long nl  = gcid->len/sizeof(int);

    coord ips[nl];
    coord uip[nl];
    coord nsg[nl];

    foreach() {
        dg0[] = 0;
        neg[] = nodata;
        if (gid[] > -1) {
            fragment interFrag;
            coord fc, gc = {x,y,z};
            PointIBM bioff;

            coord bi, ip, n = {0,0,0};
            if (cs.ibm.bi) 
                bi = cs.ibm.bi(gc);
            else {
                closest_interface (point, midPoints, cs, normals, ibalphas, &interFrag, &fc, &bioff, n);
                coord nplic = {normals.x[bioff.i, bioff.j, bioff.k], 
                               normals.y[bioff.i, bioff.j, bioff.k],
                               normals.z[bioff.i, bioff.j, bioff.k]};
                bi = boundary_int (point, nplic, nplic, interFrag.alpha, fc, cs);
                n = interpolate_normal (point, bi, bioff, midPoints, normals);
                //n = nplic;

#if 0
                for (int i = 0; i < 10; i++) {
                    closest_interface (point, midPoints, cs, normals, ibalphas, &interFrag, &fc, &bioff, n);
                    nplic = (coord){normals.x[bioff.i, bioff.j, bioff.k], 
                                    normals.y[bioff.i, bioff.j, bioff.k],
                                    normals.z[bioff.i, bioff.j, bioff.k]};
                    bi = boundary_int (point, n, nplic, interFrag.alpha, fc, cs);
                    n = interpolate_normal_lsq (point, bi, bioff, midPoints, normals);
                }
#endif
            }
            ip = image_point (bi, gc, n);
            #if 1

            coord ne = cs.ibm.normal(gc);

            double nerror = acos(clamp(dot_product(n, ne) / 
                                (magnitude_coord(n) * magnitude_coord(ne)), -1, 1));

            neg[] = nerror;
            //if (fabs(nerror) > 1e-4 && fabs(nerror) < 1e-3) 
            //if (cs.ibm.normal && !is_deep_ghost_cell(point))
                n = ne;

            normalize_sum(&n);
            #else

            if (cs.ibm.normal) {
            coord ne = cs.ibm.normal(gc);

            double nerror = acos(clamp(dot_product(n, ne) / 
                                (magnitude_coord(n) * magnitude_coord(ne)), -1, 1));
            neg[] = nerror;

            double thetae = 1e-2;
            
            normalize(&ne);

            if (false && magnitude_coord(gc) > 0.375) 
            {
              n.x = cos(thetae)*ne.x - sin(thetae)*ne.y;
              n.y = sin(thetae)*ne.x + cos(thetae)*ne.y;
            }

            normalize_sum(&n);
            }
            #endif

            dg0[] = distance_coord(bi, gc)/Delta;

            foreach_dimension()
                bis.x[] = bi.x;

            nsg[(int)gid[]] = n;
            ips[(int)gid[]] = ip;

        }
    }

    /**
    When using MPI, an image point may lay in another rank, requiring extra communication. 
    We identify these problematic points and store its index in the bips (bad image points) array.
    We also have to store the rank the IP belongs to in an array, bpid. */

#if _MPI
    Array * bips  = array_new();
    Array * bpid  = array_new();
#endif

    for (int i = 0; i < nl; i++) {

        int rank = -1;
        uip[i] = interpolate_image_point(ips[i], nsg[i], &rank);

#if _MPI
        if (uip[i].x == nodata) {
            assert (rank >= 0 && rank < npe());
            array_append(bips, &i, sizeof(int));
            array_append(bpid, &rank, sizeof(int));
        }
#else
        if (uip[i].x == nodata) {
            fprintf(stderr, "ERROR: ips[%d]=(%g, %g) nsg[%d]=(%g, %g)\n",
                i, ips[i].x, ips[i].y, i, nsg[i].x, nsg[i].y);
        }
        assert(uip[i].x != nodata);
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
    types, like MPIcoord, are not natively supported in MPI.

    TODO: can be optimized to avoid having repeat arrays. */

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
    Note we do not have to sort this array, since recv_ips is already sorted.

    TODO: use a different data type other than MPIcoord since we don't use nx, ny, and nz. */

    MPIcoord * uip_send = malloc(recv_total * sizeof(MPIcoord));

    for (int i = 0; i < recv_total; i++) {

        MPIcoord data = recv_ips[i];
        coord ip = {data.x, data.y, data.z};
        coord n = {data.nx, data.ny, data.nz};
        
        int rank = -1;
        coord uip = interpolate_image_point(ip, n, &rank);

        assert(rank == pid() && uip.x != nodata);

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

    if (uip_send) free(uip_send);
    if (uip_recv) free(uip_recv);
    if (send_ips) free(send_ips);
    if (recv_ips) free(recv_ips);
    
    array_free(bips);
    array_free(bpid);

#endif // _MPI

    array_free(gcid);

    u.x.mp = bis; u.y.mp = bis; 
#if dimension == 3
    u.z.mp = bis;
#endif

    foreach() {
        //if (is_ghost_cell(point, cs)) {
        if (gid[] > -1) {

#if 0
            if (gidd[] > 0) {
                foreach_dimension()
                    u.x[] = 0;
                continue;
            }
                
#endif

#if 0
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

            coord imageVelocity = image_velocity (point, u, ip, n, midPoints, normals, ibalphas);
#else
            coord bi = {0};
            foreach_dimension()
                bi.x = bis.x[];

            coord gc = {x, y, z};

            coord imageVelocity = uip[(int)gid[]];
            coord n = nsg[(int)gid[]];
#endif

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
                                delta += sq(2*(gc.x - bi.x));
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
    boundary((scalar *){u, bis});
}

event end_timestep(i++)
{
    fill_interface_data(); 
    update_gc_velocity();
}

#if 0
void ibm_force_distribution(vector Fmu)
{
    /**
    First count all ghost cells. */
    int gcount = 0;

    Array * gcid  = array_new();

    foreach(serial, reduction(+:gcount)) {
        gid[] = gbid[] = -1;
        if (is_ghost_cell(point, cs)) {
            gid[] = gcount;
            array_append(gcid, &gcount, sizeof(int));
            gcount++;
        }
    }

    /**
    Allocate arrays for storing the, IP, IP velocity, and normals for each GC */

    long nl  = gcid->len/sizeof(int);

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

            foreach_dimension()
                bis.x[] = bi.x;

            nsg[(int)gid[]] = n;
            ips[(int)gid[]] = ip;

        }
    }

    /**
    When using MPI, an image point may lay in another rank, requiring extra communication. 
    We identify these problematic points and store its index in the bips (bad image points) array.
    We also have to store the rank the IP belongs to in an array, bpid. */

#if _MPI
    Array * bips  = array_new();
    Array * bpid  = array_new();
#endif

    for (int i = 0; i < nl; i++) {

        int rank = -1;
        uip[i] = interpolate_image_point(ips[i], nsg[i], &rank);

#if _MPI
        if (uip[i].x == nodata) {
            assert (rank >= 0 && rank < npe());
            array_append(bips, &i, sizeof(int));
            array_append(bpid, &rank, sizeof(int));
        }
#else
        assert(uip[i].x != nodata);
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
    types, like MPIcoord, are not natively supported in MPI.

    TODO: can be optimized to avoid having repeat arrays. */

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
    Note we do not have to sort this array, since recv_ips is already sorted.

    TODO: use a different data type other than MPIcoord since we don't use nx, ny, and nz. */

    MPIcoord * uip_send = malloc(recv_total * sizeof(MPIcoord));

    for (int i = 0; i < recv_total; i++) {

        MPIcoord data = recv_ips[i];
        coord ip = {data.x, data.y, data.z};
        coord n = {data.nx, data.ny, data.nz};
        
        int rank = -1;
        coord uip = interpolate_image_point(ip, n, &rank);

        assert(rank == pid() && uip.x != nodata);

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

    if (uip_send) free(uip_send);
    if (uip_recv) free(uip_recv);
    if (send_ips) free(send_ips);
    if (recv_ips) free(recv_ips);
    
    array_free(bips);
    array_free(bpid);

#endif // _MPI

    array_free(gcid);

    foreach() {
        foreach_dimension()
            Fmu.x[] = 0
        if (cs[] > 0 && cs[] < 1) {
            coord dudn = {0,0,0};
            foreach_neighbor(1) {
                if (gid[] > -1) {
                    if (fabs(bi.x - x) < 0.5 && fabs(bi.y - y) < 0.5 && fabs(bi.z - z) < 0.5) {

                        // TODO: what if cell has more than one BI? (avg)
                        coord bi = {0,0,0};
                        foreach_dimension()
                            bi.x = bis.x[];

                        coord uip0 = uip[(int)gid[]];
                        coord ugc0 = (coord){u.x[], u.y[], u.z[]};

                        coord n = nsg[(int)gid[]];
                        coord gc = {x, y, z}, ip = image_point (bi, gc, n);
                        
                        foreach_dimension()
                            dudn.x = (uip0.x - ugc0.x)/(ip.x - gc.x);
                    }
                }
            }
            if (!empty_coord(dudn)) {
                Fmu.x[] = 
            }
        }
    }

    foreach() {
        if (is_ghost_cell(point, cs)) {

            coord gc = {x, y, z};

            coord bi = {0};
            foreach_dimension()
                bi.x = bis.x[];

            coord imageVelocity = uip[(int)gid[]];
            coord n = nsg[(int)gid[]];

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
                                delta += sq(2*(gc.x - bi.x));
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
    boundary((scalar *){u, bis});
}
#endif
