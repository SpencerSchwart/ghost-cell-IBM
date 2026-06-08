/**
Temperature module for the ghost-cell immersed boundary method. */

#include "ibm/src/my-bcg.h"
#include "ibm/src/my-diffusion.h"

scalar T[];

scalar * ibm_tracers = {T};

attribute {
    face vector diffusion;
}

event tracer_advection (i++)
{
    advection (ibm_tracers, uf, dt);
}

double image_point_value (scalar s, coord ip, coord nsg, int * rank = NULL) 
{
    double sip = nodata;
#if _OPENMP
    foreach_image_point (ip.x, ip.y, ip.z, rank, serial)
#else
    foreach_image_point (ip.x, ip.y, ip.z, rank)
#endif
        sip = interpolate_image_point_scalar (point, s, ip, nsg, midPoints, normals, ibalphas);

    return sip;
}

void update_ghost_values()
{
    /**
    MPI boundaries MUST be updated before interpolation. */
#if _MPI
    mpi_boundary_update (ibm_tracers);
#endif

    /**
    First count all ghost cells. */
    int gcount = 0;

    Array * gcid  = array_new();

    foreach(serial, reduction(+:gcount)) {
        gid[] = -1;
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
    coord nsg[nl];
    double sips[nl];

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
                n = interpolate_normal_lsq (point, bi, bioff, midPoints, normals);
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
        sips[i] = image_point_value(T, ips[i], nsg[i], &rank);
        if (sips[i] == nodata) {
            fprintf(stderr, "ERROR: ips[%d]=(%g, %g) nsg[%d]=(%g, %g)\n",
                i, ips[i].x, ips[i].y, i, nsg[i].x, nsg[i].y);
        }
        assert(sips[i] != nodata);
    }

    array_free(gcid);

    foreach() {
        if (gid[] > -1) {
            coord bi = {0};
            foreach_dimension()
                bi.x = bis.x[];

            coord gc = {x, y, z};

            double sval = sips[(int)gid[]];

            bool bctype[2] = {false, false};
            double vb = T.boundary[immersed] (point, point, T, bctype);
            bool dirichlet = bctype[0], navier = bctype[1];
            if (dirichlet) {
                if (navier) {
                    double delta = 0;
                    foreach_dimension()
                        delta += sq(2*(gc.x - bi.x));
                    delta = sqrt(delta)/(2*Delta);
                    vb /= Delta;

                    T[] = -sval*(delta - vb)/(delta + vb);
                }
                else {
                    T[] = 2*vb - sval;
                }
            }
            else { // neumann
                    T[] = sval;
            }
       }
       else if (cs[] == 0) {
           T[] = 0.;
       }
    }
    boundary(ibm_tracers);
}


event tracer_diffusion (i++)
{
    diffusion(T, dt, T.diffusion);

    fill_interface_data(); 
    update_ghost_values();
}

double nusselt_number(scalar T, scalar nu)
{
    double avgnu = 0;
    double tarea = 0;

    foreach(reduction(+:avgnu) reduction(+:tarea)) {
        nu[] = nodata;
        if (cs[] > 0 && cs[] < 1) {
            coord n = {normals.x[], normals.y[]}, pint, ip;

            double area = pow(Delta, dimension - 1)*plane_area_center(n, ibalphas[], &pint);
            
            pint = local_to_global_coord(point, pint);

            normalize2(&n);
            foreach_dimension()
                ip.x = pint.x + Delta*n.x;
             
            double Tip = interpolate_image_point_scalar (point, T, ip, n, midPoints, normals, ibalphas);
            
            bool dirichlet = false;
            double Tb = T.boundary[immersed](point, point, T, &dirichlet);
            double dTdn = dirichlet? (Tip - Tb)/Delta: Tb;
            
            nu[] = -dTdn;

            avgnu += dTdn*area;
            tarea += area;
        }
    }

    return tarea? -avgnu/tarea: 0;
}

