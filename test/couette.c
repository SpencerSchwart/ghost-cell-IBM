/**
# Couette flow between rotating cylinders

We test immersed boundaries by solving the (Stokes) Couette flow
between two rotating cylinders. */

#define GCV 1.09e-4

#include "grid/quadtree.h"
#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "view.h"

#include "navier-stokes/perfs.h"

int bc_type = 0;           // boundary condition type (dirichlet = 0, neumann = 1, navier = 2)
double lambda = 0.2;       // slip length for navier boundary condition

double ri = 0.25, ro = 0.5; // radius of inner and outter cylinder

bool append = false;        // append to existing outfile.

coord solid_normal (coord p)
{
    coord n = {p.x, p.y};
    double d = sqrt(sq(p.x) + sq(p.y));

    if (d > 0.5*(ri + ro))
        foreach_dimension()
            n.x *= -1;

    normalize_sum(&n);
    return n;
}

coord solid_bi (coord p)
{
    double d = sqrt(sq(p.x) + sq(p.y));
    double R = d > 0.5*(0.25 + 0.25)? 0.5: 0.25;
    
    return (coord){R*p.x/d, R*p.y/d};
}

int main(int argc, char* argv[])
{
  mu = fm;

  size (1.25);
  DT = 1e-2;

  origin (-L0/2., -L0/2.);
  stokes = true;

  //cs.ibm.bi = solid_bi;
  cs.ibm.normal = solid_normal;

  if (argc > 1)
    bc_type = atoi (argv[1]);
  if (argc > 2)
    lambda = atof (argv[2]);
  if (argc > 3) {
    N = atoi (argv[3]);
    append = true;
    run();
    return 0;
  }
 
  if (bc_type == 4) {
    for (bc_type = 0; bc_type < 3; bc_type++) {
      for (N = 16; N <= 128; N *= 2) {
        run();
      }
    }
  }
  else {
    for (N = 16; N <= 512; N *= 2) 
      run();
  }
}

scalar un[];

double powerlaw_noslip(double r)
{
  return (r*(sq(ro/r) - 1.)/(sq(ro/ri) - 1.));
}

double powerlaw_slip(double r)
{
  return sq(ri)/(sq(ri) + sq(ro))*(sq(ro)/r + r);
}

double powerlaw_navier(double r)
{
  return (r*(ro-lambda) + sq(ro)*(-lambda-ro)/r)*(sq(ri)/(sq(ri)*(ro-lambda)+sq(ro)*(-lambda-ro)));
}

double (*theory) (double);

FILE * fout;

event init (t = 0) {
  
  /**
  The geometry of the embedded boundary is defined as two cylinders of
  radii 0.5 and 0.25. */

  solid (cs, fs, difference (sq(ro) - sq(x) - sq(y),
			     sq(ri) - sq(x) - sq(y)));

  //event("update_metric");

  /**
  The outer cylinder is fixed and the inner cylinder is rotating with
  an angular velocity unity. */
  
  /**
  We initialize the reference velocity field. */
  
  foreach()
    un[] = u.y[];


  /**
  Initalize the wall boundary conditions. Each case has a no penetration condition, 
  i.e. un = 0. */

  u.n[immersed] = dirichlet(0);

  char nameo[80];
  switch (bc_type) {
    case 0:
      u.t[immersed] = sqrt(x*x + y*y) > 0.5*(ro + ri)? dirichlet(0): dirichlet(distance(x, y));
      sprintf(nameo, "outs/%g-noslip-out-%d", GCV, N);
      theory = powerlaw_noslip;
      break;
    case 1:
      u.t[immersed] = sqrt(x*x + y*y) > 0.5*(ro + ri)? neumann(0): dirichlet(distance(x, y));
      theory = powerlaw_slip;
      sprintf(nameo, "outs/%g-slip-out-%d", GCV, N);
      break;
    case 2:
      u.t[immersed] = sqrt(x*x + y*y) > 0.5*(ro + ri)? navier(lambda): dirichlet(distance(x, y));
      theory = powerlaw_navier;
      sprintf(nameo, "outs/%g-navier-%g-out-%d", GCV, lambda, N);
      break;
  }

  fout = fopen(nameo, "w");
}

/**
We look for a stationary solution. */

scalar e0[], e[], egc[], uthetat[], uthetas[];

event logfile (i++; i <= 10000) {
 foreach() {
    e0[] = e[];
    e[] = egc[] = uthetas[] = uthetat[] = nodata;
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    if (cs[] > GCV) {
      uthetas[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      uthetat[] = theory(r);
      e[]  = fabs(fabs(uthetas[]) - fabs(theory(r)));
    }
    else if (is_ghost_cell(point, cs)) {
      uthetas[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      uthetat[] = theory(r);
      egc[]  = fabs(fabs(uthetas[]) - fabs(theory(r)));
    }
  }

  norm2 en   = normf2(e);
  norm2 engc = normf2(egc);

  norm2 nnorm = normf2(neg);
  stats nstat = statsf(neg);

  double du = change (u.y, un);

  fprintf(fout, "%d %g %g %g %g %g %g %g %g %g\n", i, t, du, nnorm.avg, nnorm.max, nstat.stddev, engc.rms, engc.max, en.rms, en.max);
  fflush(fout);
 
  /**
  For some reason, some cases using openMP cannot get below du_max < 1e-8 */
#if _OPENMP
  if (i > 0 && du < 1e-6)
#else
  if (i > 0 && du < 1e-8)
#endif
    return 1; /* stop */
}

/**
We compute error norms and display the angular velocity, pressure and
error fields using bview. */

event profile (t = end)
{
  scalar utheta[], efl[], edgc[], engc[];
  foreach() {
    double theta = atan2(y, x), r = sqrt(x*x + y*y);
    e[] = efl[] = egc[] = edgc[] = engc[] = nodata;
    if (cs[] > GCV) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      efl[] = fabs(utheta[] - theory(r));

      if (cs[] > 0.5) 
        e[]  = fabs(utheta[] - theory(r));
    }
    else if (is_ghost_cell(point, cs)) {
      utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
      egc[]  = fabs(utheta[] - theory(r));
      
      bool inside = true;
      foreach_direct_neighbor()
        if (gc[])
            inside = false;

      if (inside) 
          edgc[] = fabs(utheta[] - theory(r));
      else 
          engc[] = fabs(utheta[] - theory(r));
    }
    else
      p[] = utheta[] = nodata;
  }

  char name[80], namei1[80], namei2[80], namei3[80], namei4[80], namei5[80], namei6[80], namei7[80], namei8[80], namei9[80];

  switch (bc_type) {
    case 0:
      sprintf(name, "data/%g-noslip.dat", GCV);
      sprintf(namei1, "imgs/%g-u-noslip-%d.png", GCV, N);
      sprintf(namei2, "imgs/%g-p-noslip-%d.png", GCV, N);
      sprintf(namei3, "imgs/%g-e-noslip-%d.png", GCV, N);
      sprintf(namei4, "imgs/%g-d-noslip-%d.png", GCV, N);
      sprintf(namei5, "imgs/%g-un-noslip-%d.png", GCV, N);
      sprintf(namei6, "imgs/%g-eu-noslip-%d.png", GCV, N);
      sprintf(namei7, "imgs/%g-en0-noslip-%d.png", GCV, N);
      sprintf(namei8, "imgs/%g-en-noslip-%d.png", GCV, N);
      sprintf(namei9, "imgs/%g-cell-noslip-%d.png", GCV, N);
      break;
    case 1:
      sprintf(name, "data/%g-slip.dat", GCV);
      sprintf(namei1, "imgs/%g-u-slip-%d.png", GCV, N);
      sprintf(namei2, "imgs/%g-p-slip-%d.png", GCV, N);
      sprintf(namei3, "imgs/%g-e-slip-%d.png", GCV, N);
      sprintf(namei4, "imgs/%g-d-slip-%d.png", GCV, N);
      sprintf(namei5, "imgs/%g-nu-slip-%d.png", GCV, N);
      sprintf(namei6, "imgs/%g-eu-slip-%d.png", GCV, N);
      sprintf(namei7, "imgs/%g-en0-slip-%d.png", GCV, N);
      sprintf(namei8, "imgs/%g-en-slip-%d.png", GCV, N);
      sprintf(namei9, "imgs/%g-slip-slip-%d.png", GCV, N);
      break;
    case 2:
      sprintf(name, "data/%g-navier-%g.dat", GCV, lambda);
      sprintf(namei1, "imgs/%g-u-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei2, "imgs/%g-p-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei3, "imgs/%g-e-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei4, "imgs/%g-d-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei5, "imgs/%g-nu-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei6, "imgs/%g-eu-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei7, "imgs/%g-en0-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei8, "imgs/%g-en-navier-%g-%d.png", GCV, lambda, N);
      sprintf(namei9, "imgs/%g-cell-navier-%g-%d.png", GCV, lambda, N);
      break;
  }

  static FILE * fp = append? fopen(name, "a"): fopen(name, "w");

  /**
  Calculate the divergence in 1. all cells 2. fluid cells, 3. (solid) interfacial cells 
  4. small cells (solid volume > 0.5). Also interpolate to find the velocity error on the interface,
  where us.x = us.n and us.y = us.t. */

  scalar div[], divf[], divsi[], divsc[];
  vector us[];
  scalar ue[];          // tangential velocity error on solid
  scalar dn0[], dn1[];  // distance between GC to IP (not normalize and normalized by cell size)
  scalar nse[], nse0[]; // normal error
  scalar dse[];         // distance from BI to midpoint

  foreach() {
    dse[] = nse0[] = nse[] = dn0[] = dn1[] = ue[] = div[] = divf[] = divsi[] = divsc[] = nodata;

    foreach_dimension()
        us.x[] = nodata;

    if (cs[]) {

      foreach_dimension()
        div[] += fs.x[1]*uf.x[1] - fs.x[]*uf.x[];

      div[] *= Delta;

      if (cs[] >= 1) { // pure fluid cell
        divf[] = 0;
        foreach_dimension()
          divf[] += fs.x[1]*uf.x[1] - fs.x[]*uf.x[];
        divf[] *= Delta;
      }

      if (cs[] < 1 && cs[] > 0) { // sold interface cells
        divsi[] = 0;
        foreach_dimension()
          divsi[] += fs.x[1]*uf.x[1] - fs.x[]*uf.x[];
          
        divsi[] *= Delta;

        coord n = facet_normal (point, cs, fs), pint;
        coord ts = {-n.y, n.x};
        double alpha = plane_alpha (cs[], n);
        plane_area_center (n, alpha, &pint);
        pint = local_to_global_coord(point, pint);

        coord ui;
        ui.x = interpolate_linear (point, u.x, pint.x, pint.y, 0);
        ui.y = interpolate_linear (point, u.y, pint.x, pint.y, 0);

        normalize(&n); normalize(&ts);

        us.x[] = dot_product(n, ui);
        us.y[] = dot_product(ts, ui);

        double r = sqrt(sq(pint.x) + sq(pint.y));

        if (r > 0.5*(ro + ri))
          ue[] = fabs(fabs(us.y[]) - fabs(theory(r)));

      }
      if (cs[] < 0.5 && cs[] > 0) { // small cells
        divsc[] = 0;
        foreach_dimension()
          divsc[] += fs.x[1]*uf.x[1] - fs.x[]*uf.x[];
        divsc[] *= Delta;
      }
    }

      // calculate distance from GC to IP
      if (is_ghost_cell(point, cs)) {
        fragment interFrag;
        coord fluidCell, ghostCell = {x,y,z}, n = {0,0,0};
        PointIBM bioff;

        coord tmp = closest_interface (point, midPoints, cs, normals, ibalphas, &interFrag, &fluidCell, &bioff, n);
        coord boundaryInt = boundary_int (point, interFrag.n, interFrag.n, interFrag.alpha, fluidCell, cs);
        coord imagePoint = image_point (boundaryInt, ghostCell);
 
        dn0[] = distance_coord(imagePoint, ghostCell);
        dn1[] = distance_coord(imagePoint, ghostCell)/Delta;

        coord m0 = interFrag.n;
        coord m  = interpolate_normal (point, boundaryInt, bioff, midPoints, normals);
        coord me = solid_normal(ghostCell);

        nse0[] = acos(clamp(dot_product(m0, me) / 
                     (magnitude_coord(m0) * magnitude_coord(me)), -1, 1));
        nse[] = acos(clamp(dot_product(m, me) / 
                    (magnitude_coord(m) * magnitude_coord(me)), -1, 1));
        dse[]  = distance_coord(tmp, boundaryInt);
      }
  }

  norm2 n    = normf2(e);     // fluid cell error (cs > 0.5)
  norm2 ngc  = normf2(egc);   // all ghost cell error
  norm2 ndgc = normf2(edgc);  // deep ghost cell error
  norm2 nngc = normf2(engc);  // normal ghost cell error
  norm2 nfl  = normf2(efl);   // fluid cell error (cs > GCV)
  norm2 ma   = normf2(div);   // mass error
  norm2 mf   = normf2(divf);  // mass error (only pure fluid cells)
  norm2 ms   = normf2(divsi); // mass error (only interface cells)
  norm2 msc  = normf2(divsc); // mass error (only small cells)

  norm2 nun  = normf2(us.x);  // normal velocity on wall (should be 0 from no penetration BC)
  norm2 nut  = normf2(ue);    // tangential velocity error

  norm2 nd0   = normf2(dn0);    // distance from image point
  norm2 nd1   = normf2(dn1);    // distance from image point (normalized)

  norm2 nnse0 = normf2(nse0); 
  norm2 nnse  = normf2(nse);
  norm2 ndse  = normf2(dse);

  norm2 nnorm = normf2(neg);
  stats nstat = statsf(neg);

  if (N == 16 && !append)
    fprintf(fp, "0:N 1:n.avg 2:n.rnms 3:n.max 4:n.l2 5:i 6:gp.i 7:mgp.nrelax 8:mgu.i 9:mgu.nrelax 10:lambda 11:ma.avg 12:ma.rms 13:ma.max 14:mf.avg" 
                " 15:mf.rms 16:mf.max 17:ms.avg 18:ms.rms 19:ms.max 20:msc.avg 21:msc.rms 22:msc.max 23:nun.avg 24:nun.rms 25:nun.max 26:nun.l2" 
                " 27:nut.avg 28:nut.rms 29:nut.max 30:nut.l2 31:nd0.avg 32:nd0.rms 33:nd0.max 34:nd0.l2 35:nd1.avg 36:nd1.rms 37:nd1.max 38:nd1.l2"
                " 39:ngc.rms 40:ngc.max 41:ngc.l2 42:ndgc.rms 43:ndgc.max 44:ndgc.l2 45:nngc.rms 46:nngc.max 47:nngc.l2 48:nfl.rms 49:nfl.max 50:nfl.l2"
                " 51:nnse0.avg 52:nnse0.rms 53:nnse0.max 54:nnse.avg 55:nnse.rms 56:nnse.max 57:ndse.avg 58:ndse.rms 59:ndse.max"
                " 60:nnorm.avg 61:nnorm.rms 62:nnorm.max 63:nstat.stddev\n");

  fprintf (fp, "%d %.7g %.7g %.7g %.7g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   N, n.avg, n.rms, n.max, n.l2, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, lambda,
       ma.avg, ma.rms, ma.max, mf.avg, mf.rms, mf.max, ms.avg, ms.rms, ms.max,
       msc.avg, msc.rms, msc.max, nun.avg, nun.rms, nun.max, nun.l2, nut.avg, nut.rms, nut.max, nut.l2,
       nd0.avg, nd0.rms, nd0.max, nd0.l2, nd1.avg, nd1.rms, nd1.max, nd1.l2, ngc.rms, ngc.max, ngc.l2,
       ndgc.rms, ndgc.max, ndgc.l2, nngc.rms, nngc.max, nngc.l2, nfl.rms, nfl.max, nfl.l2,
       nnse0.avg, nnse0.rms, nnse0.max, nnse.avg, nnse.rms, nnse.max, ndse.avg, ndse.rms, ndse.max,
       nnorm.avg, nnorm.rms, nnorm.max, nstat.stddev);

  fflush(fp);
  if (N == 512)
    fclose(fp);

  /**
  Output images of the velocity, pressure, and error fields */
  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("utheta", spread = -1);
  save (namei1);

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("p", spread = -1);
  save (namei2);

  draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
  squares ("e", min = -n.max, max = n.max, cbar=true);
  save (namei3);

  draw_vof ("cs", "fs");
  squares ("dn1", min = 1, max = 3, cbar=true);
  save (namei4);

  draw_vof ("cs", "fs");
  squares ("us.x", min = -nun.max, max = nun.max, cbar=true);
  save (namei5);

  draw_vof ("cs", "fs");
  squares ("ue", min = -nut.max, max = nut.max, cbar=true);
  save (namei6);

  draw_vof ("cs", "fs");
  squares ("nse0", min = -nnse0.max, max = nnse0.max, cbar=true);
  save (namei7);

  draw_vof ("cs", "fs");
  squares ("nse", min = -nnse.max, max = nnse.max, cbar=true);
  save (namei8);

  view(width=1600, height=1600);
  draw_vof ("cs", "fs", filled=-1, lw = 2, lc={0.5, 0.5, 0.5}, fc={0.796875, 0.796875, 0.796875});
  cells(lw = 0.5);
  save (namei9);

  /**
  Print velocity profiles for N = 32 run. */
  if (N == 32) {
    switch (bc_type) {
      case 0:
        sprintf(name, "profiles/%g-noslip-prof.dat", GCV);
        break;
      case 1:
        sprintf(name, "profiles/%g-slip-prof.dat", GCV);
        break;
      case 2:
        sprintf(name, "profiles/%g-navier-prof-%g.dat", GCV, lambda);
        break;
  }

    FILE * fpv = fopen(name, "w");
    foreach() {
      double theta = atan2(y, x), r = sqrt(x*x + y*y);
      if (fabs(utheta[]) < nodata)
          fprintf (fpv, "%g %g %g %g %g %g %g\n",
	           r, theta, u.x[], u.y[], p[], utheta[], e[]);
    }
    fclose(fpv);
  }
  
  /**
  Print boundary intercepts for different grids */
  sprintf(name, "data/%g-bi-points-%d", GCV, N);
  FILE * fpb = fopen(name, "w");
  fprintf(fpb, "0:x 1:y 2:bi.x 3:bi.y 4:ip.x 5:ip.y 6:bie.x 7:bie.y 8:ipe.x 9:ipe.y 10:m.x 11:m.y 12:me.x 13:me.y 14:nerror 15:dgcip 16:deep 17:egc\n");
  foreach() {
    if (is_ghost_cell(point, cs)) {
      fragment interFrag;
      coord fluidCell, gc = {x, y}, n = {0,0};
      PointIBM bioff;

      closest_interface (point, midPoints, cs, normals, ibalphas, &interFrag, &fluidCell, &bioff, n);
      coord bi = boundary_int (point, interFrag.n, interFrag.n, interFrag.alpha, fluidCell, cs);
      coord ip = image_point (bi, gc);

      // Calculate exact BI and IP
      double d = sqrt(sq(x) + sq(y));
      double R = d > 0.5*(ri + ro)? ro: ri;

      coord bie = {R*x/d, R*y/d};
      coord ipe = image_point(bie, gc);

      coord m  = interpolate_normal (point, bi, bioff, midPoints, normals);
      coord me = solid_normal(gc);
      double nerror = acos(clamp(dot_product(m, me) / 
                          (magnitude_coord(m) * magnitude_coord(me)), -1, 1));

      double dgcip = distance_coord(ip, gc)/Delta;

      bool deep = is_deep_ghost_cell(point);

      ip = image_point(bi, gc, m);

      fprintf(fpb, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g\n", 
        x, y, bi.x, bi.y, ip.x, ip.y, bie.x, bie.y, ipe.x, ipe.y, m.x, m.y, me.x, me.y, nerror, dgcip, deep, egc[]);
    }
  }
  fclose(fpb);

  sprintf(name, "interfaces/interface-%d", N);
  fpb = fopen(name, "w");
  output_facets(cs, fpb, fs);
}

