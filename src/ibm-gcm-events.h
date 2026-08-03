
vector bis[];

scalar ggid[];

void update_metric();

event init(i = 0)
{
    ggid.restriction = ghost_cell_id_restriction;
    ggid.refine = ghost_cell_id_restriction;
    update_metric();

    boundary((scalar *){cs});
}

void fill_interface_data (scalar s, scalar alphas, vector ns, vector mps) {
    foreach() {
        coord mp, n;
        if (s[] > 0 && s[] < 1) {
            centroid_point (point, s, fs, &mp, &n, &alphas[]);
            foreach_dimension() {
                mps.x[] = mp.x;
                ns.x[] = -n.x;
            }
            alphas[] *= -1;
        }
        else if (s[] == 1 && empty_neighbor (point, &mp, &n, s)) {
            foreach_dimension() {
                mps.x[] = mp.x;
                ns.x[] = -n.x;
            }
            alphas[] = 0;
        }
        else {
            foreach_dimension() {
                mps.x[] = 0;
                ns.x[] = 0;
            }
            alphas[] = 0;
        }
    }

    boundary({mps, ns, alphas});
}

#if _MPI
typedef struct {
    int id;    // ID in array
    double x, y, z;
    double nx, ny, nz;
} MPIcoord;

#endif

// TODO: doesn't work in openmp?
coord interpolate_image_point (coord ip, coord nsg, vector mps, vector ns, scalar alphas, int * rank = NULL) 
{
    coord uip = {nodata};
#if _OPENMP
    foreach_image_point (ip.x, ip.y, ip.z, rank, serial)
#else
    foreach_image_point (ip.x, ip.y, ip.z, rank)
#endif
        uip = interpolate_vector_ibm (point, u, ip, nsg, mps, ns, alphas);

    return uip;
}


void update_gc_velocity (scalar cs, scalar alphas, vector ns, vector mps)
{
    /**
    MPI boundaries MUST be updated before interpolation. */
    boundary((scalar *){u});

    /**
    First count all ghost cells. gid index's the ghost cells starting at 1, not
    0. This is because for plane ghost cells, gid = 0. */
    scalar gid[];
    int gcount = 0;

    Array * gcid  = array_new();

    foreach(serial, reduction(+:gcount)) {
        ggid[] = gid[] = -1;
        if (is_ghost_cell(point, cs)) {
            gcount++;
            ggid[] = gid[] = gcount;
            array_append(gcid, &gcount, sizeof(int));
        }
    }

    /**
    Allocate arrays for storing the, IP, IP velocity, and normals for each GC */

    long nl  = gcid->len/sizeof(int);

    coord ips[nl];
    coord uips[nl];
    coord nsg[nl];

    foreach() {
        
        foreach_dimension()
            bis.x[] = 0;

        if (gid[] > 0) {
            coord fc, gc = {x,y,z};
            PointIBM bioff;

            coord bi, ip, n = {0,0,0};
            if (cs.ibm.bi) 
                bi = cs.ibm.bi(gc);
            else {
                double talpha = 0;
                closest_interface (point, mps, cs, ns, alphas, &talpha, &fc, &bioff, n);
                coord nplic = {ns.x[bioff.i, bioff.j, bioff.k], 
                               ns.y[bioff.i, bioff.j, bioff.k],
                               ns.z[bioff.i, bioff.j, bioff.k]};
                bi = boundary_int (point, nplic, nplic, talpha, fc, cs);
                n = interpolate_normal_lsq (point, bi, bioff, mps, ns);
            }
            if (cs.ibm.normal)
               n = cs.ibm.normal(gc);

            ip = image_point (bi, gc, n);

            foreach_dimension()
                bis.x[] = bi.x;

            nsg[(int)gid[]-1] = n;
            ips[(int)gid[]-1] = ip;

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
        uips[i] = interpolate_image_point(ips[i], nsg[i], mps, ns, alphas, &rank);

        if (rank < 0)   // image point is in boundary
            continue;
#if _MPI
        if (uips[i].x == nodata) {
            assert (rank >= 0 && rank < npe());
            array_append(bips, &i, sizeof(int));
            array_append(bpid, &rank, sizeof(int));
        }
#else
        if (uips[i].x == nodata) {
            fprintf(stderr, "ERROR: ips[%d]=(%g, %g) nsg[%d]=(%g, %g)\n",
                i, ips[i].x, ips[i].y, i, nsg[i].x, nsg[i].y);
        }
        assert(uips[i].x != nodata);
#endif
    }

    /**
    If "bad" image points exists, we count the total number for each rank and store it
    in an array whose indices correspond to all active processes. We then communicate this
    array to processors. */
#if 1
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
        coord uip = interpolate_image_point(ip, n, mps, ns, alphas, &rank);

        assert(rank >= 0);

        /**
        In rare cases, the initial rank determination of the image point may be
        from a parent cell who's children belong to another process, partly caused
        from each tree having a different maximum level (depth()). This requires
        another round of communication, potentially more(?), but it is not 
        implemented yet. Instead, we rely on the image point velocity being 
        calculated later on within the ghost cell. */
        //if (rank != pid())

        //assert(rank == pid());
        //assert(uip.x != nodata);

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
        uips[data.id] = (coord){data.x, data.y, data.z};
    }

    if (uip_send) free(uip_send);
    if (uip_recv) free(uip_recv);
    if (send_ips) free(send_ips);
    if (recv_ips) free(recv_ips);
    
    array_free(bips);
    array_free(bpid);

#endif // _MPI
#endif
    array_free(gcid);

    u.x.mp = bis; u.y.mp = bis; 
#if dimension == 3
    u.z.mp = bis;
#endif

    foreach() {
        if (gid[] > 0) {
            coord bi = {0};
            foreach_dimension()
                bi.x = bis.x[];

            coord gc = {x, y, z};

            coord uip = uips[(int)gid[]-1];
            coord n   = nsg[(int)gid[]-1];

            /**
            It is possible that the image point lies in a domain boundary or it
            requires another round of MPI communication. A contingency procedure
            has not (yet) been implemented for either of these cases. Instead,
            we calculate the image point velocity within the ghost cell stencil.
            This *may* cause problems if the IP interpolation requires nodes
            outside the local 5x5 stencil, though it seems ok in testing. */

            if (uip.x == nodata) {
                coord ip = image_point (bi, gc, n);
                uip = interpolate_vector_ibm (point, u, ip, n, mps, ns, alphas);
            }

            /**
            Project the image point velocity to the local wall coordinate 
            system. */
            coord t1, t2;
            normal_and_tangents (&n, &t1, &t2);

            coord uipp = {dot_product(uip, n),
                          dot_product(uip, t1),
                          dot_product(uip, t2)};

            coord ugc = {0,0,0};

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

                        ugc.x = -uipp.x*(delta - vb)/(delta + vb);
                    }
                    else {
                        ugc.x = 2*vb - uipp.x;
                    }
                }
                else { // neumann
                        ugc.x = uipp.x;
                }
            }

            double gcn = ugc.x, gct1 = ugc.y, gct2 = ugc.z;
            foreach_dimension()
                u.x[] = gcn*n.x + gct1*t1.x + gct2*t2.x;
       }
       else if (cs[] == 0) {
           foreach_dimension() 
               u.x[] = 0.;
       }
    }
    boundary((scalar *){u, bis});
}


/**
AMR can messup velocity inside the solid boundary, which can sometimes cause 
crashes. Here, we force all solid cells that are ghost cells to have a trivial 
solution. This can probably be avoid by having a more robust refinement/coarsening
treatment. 

This function also resets the metric fields (gc and gcf) and the ghost cell ID
field (gid). */

void reset_ibm()
{
    foreach() {
        if (((!cs[] || !csp[]) && !is_ghost_cell(point, cs))) {
            ggid[] = -1;
            p[] = pf[] = 0;
            foreach_dimension()
                u.x[] = 0;
        }
        if (ggid[] > -1 && (!cs[] || !csp[]))
            p[] = pf[] = 0;

        ggid[] = -1;
    }
}


/**
centered.h initializes the face velocity field without consdering the solid 
boundary. This can cause uf != 0 on faces completely occupied by the immersed
boundary, specifically when restarting from a dump file. We call this event at
the start to avoid this problem. */

event properties (i = 0)
{
    foreach_face()
        if (!fm.x[] || !fs.x[])
            uf.x[] = 0;
}

/**
This event imposes the values in the ghost cells depending on the imposed 
boundary conditions. This takes place at the end of the events sequence,
specifically after the projection step. 

If AMR/Trees are used, we need to make sure the ghost cells aren't refined or
coarsened poorly. To avoid this, we simply wait until the grid is adapted, then
calculate and set the ghost cell values.

TODO: update_metric is called twice within a timeloop. Fix this!*/

void update_metric ()
{
    // update metrics considering immersed boundary
    foreach() {
        if (cs[] > GCV && csp[] > GCV) // fluid cell
            gc[] = 1;
        else            // ghost or solid cell
            gc[] = 0;
    }

    set_dirty_stencil (cs);
    foreach_face() {
       if ((cs[] > 0 && csp[] > 0) || (cs[-1] > 0 && csp[-1] > 0) || is_ghost_cell (point, cs))
            gcf.x[] = 1;
       else
            gcf.x[] = 0;

       if (!(csp[] && csp[-1]))
         fs.x[] = 0;
        
    }

#if AXI
    cm_update (cm, gc);
    fm_update (fm, gcf);
    boundary ({fm, cm, gcf, gc});
#else
    boundary({gcf, gc});
#endif
}


void update_symmetry(scalar s, scalar * fields)
{
    foreach() {
        if (s[] == 0 && cs[]) {
            bool is_ghost = false;
            //coord ufc = {0,0,0};
            coord fc = {0,0,0};
            foreach_direct_neighbor() 
                if (s[] == 1 && cs[]) {
                    is_ghost = true;
                    fc = (coord){x, y, z};
                    break;
                 }
            
            if (!is_ghost) {
                ggid[] = -1;
                continue;
            }

            ggid[] = 0;

            PointIBM off = {sign2(fc.x - x), sign2(fc.y - y), sign2(fc.z - z)};

            for (scalar s in fields) {
                bool bctype[2] = {false, false};
                double vb = s.boundary[cut] (point, point, s, bctype);
                bool dirichlet = bctype[0];
                if (dirichlet)
                    s[] = 2.*vb - s[off.i,off.j,off.k];
                else
                    s[] = s[off.i,off.j,off.k];
            }
        }
    }
}

event end_timestep(i++)
  update_symmetry(csp, csp.depends);

#if TREE
event adapt(i++)
{
    /**
    Reset the metric fields and anything that could've been messed up by AMR. */
    reset_ibm();

    update_metric();

    /**
    Loop through each scalar field associated with an immersed boundary and
    apply boundary conditions. */
    vector ns[], mps[];
    scalar alphas[];
    fill_interface_data(cs, alphas, ns, mps); 
    update_gc_velocity(cs, alphas, ns, mps);
    
    update_symmetry(csp, csp.depends);
}
#else
event end_timestep(i++)
{
    fill_interface_data(); 
    update_gc_velocity();
}
#endif

