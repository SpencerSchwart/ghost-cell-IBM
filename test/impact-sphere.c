/**
Droplet impacting on a sphere in axisymmetric coordinates */

#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-axi.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"

#include "view.h"
#include "tag.h"
#include "navier-stokes/perfs.h"

int maxlevel = 12;
int minlevel = 6;

double xd0 = 1.9;       // initial drop position
double ud0 = 1;         // initial drop velocity
double df = 1;          // drop diameter
double ds = 2.7;       // sphere diameter
double gravity = 9.81;

bool no_bubbles = false;

double tout = 0.05;
double tend = 12;
double trestart = 1000;
double uemax = 1e-3;

double cflg = 0.5;      // CFL number

u.t[immersed] = dirichlet(0);
u.n[immersed] = dirichlet(0);

u.n[right] = neumann(0);
p  [right] = dirichlet(0);
pf [right] = dirichlet(0);

u.n[left] = neumann(0);
p  [left] = dirichlet(0);
pf [left] = dirichlet(0);

int main(int argc, char * argv[]) {

  double l0 = 1;
  
  if (argc > 1)
    l0 = atof(argv[1]);
  if (argc > 2)
    maxlevel = atoi(argv[2]);
  if (argc > 3)
    minlevel = atoi(argv[3]);
  if (argc > 4)
    rho1 = atof(argv[4]);
  if (argc > 5)
    rho2 = atof(argv[5]);
  if (argc > 6)
    mu1 = atof(argv[6]);
  if (argc > 7)
    mu2 = atof(argv[7]);
  if (argc > 8)
    f.sigma = atof(argv[8]);
  if (argc > 9)
    gravity = atof(argv[9]);
  if (argc > 10)
    ud0 = atof(argv[10]);
  if (argc > 11)
    f.wetting.theta_s = atof(argv[11]);
  if (argc > 12)
    xd0 = atof(argv[12]);
  if (argc > 13)
    df = atof(argv[13]);
  if (argc > 14)
    ds = atof(argv[14]);
  if (argc > 15)
    tend = atof(argv[15]);
  if (argc > 16)
    tout = atof(argv[16]);
  if (argc > 17)
    trestart = atof(argv[17]);
  if (argc > 18)
    no_bubbles = (bool)atoi(argv[18]);
  if (argc > 19)
    DT = atof(argv[19]);
  if (argc > 20)
    f.wetting.dynamic = f.wetting.hysteresis_only = (bool)atoi(argv[20]);
  if (argc > 21)
    f.wetting.theta_a = atof(argv[21]);
  if (argc > 22)
    f.wetting.theta_r = atof(argv[22]);
  if (argc > 23)
    cflg = atof(argv[23]);

  TOLERANCE = 1e-4;

  double hmin = l0/(1 << maxlevel);

  size(l0);
  origin (-l0/2. -0.1*hmin, 0);
  init_grid (1 << minlevel);

  run();
}

double v0 = 0;
event init (t = 0)
{
  CFL = cflg;

  char name[80];
  sprintf(name, "dumps/dump-%g", trestart);
  if (!restore(name)) {
    scalar cs1[], f1[];
    face vector fs1[];

    /**
    Iteratively refine the initial mesh around the sphere and droplet. */

    int count = 0;
    astats st;
    do {
      solid (cs1, fs1, (sq(x) + sq(y) - sq(ds*0.5)));
      fraction (f1, - (sq(x - xd0) + sq(y) - sq(df*0.5)));
      st = adapt_wavelet ({cs1, f1}, (double[]){1e-5, 1e-5}, maxlevel, minlevel);
      count++;
    } while ((st.nc || st.nf) && count < 40);

    solid (cs, fs, (sq(x) + sq(y) - sq(ds*0.5)));
    fractions_cleanup(cs, fs, 1e-3);

    fraction (f, - (sq(x - xd0) + sq(y) - sq(df*0.5)));
    fraction (ch, - (sq(x - xd0) + sq(y) - sq(df*0.5)));

    event("update_metric");

    foreach()
      u.x[] = -ud0*f[];

    boundary(all);
  }
  else {
    solid (cs, fs, (sq(x) + sq(y) - sq(ds*0.5)));
    fractions_cleanup(cs, fs, 1e-3);
    event("update_metric");

    boundary(all);
    restriction(all);
    foreach(reduction(+:v0))
      v0 += f[]*sq(Delta)*cm[];
  }
  fprintf(stderr, "l0=%g maxl=%d minl=%d r1=%g r2=%g m1=%g m2=%g sigma=%g\n"
                  "g=%g ud0=%g theta_s=%g xd0=%g df=%g ds=%g tend=%g tout=%g trestart=%g nobub=%d DT=%g\n",
  	              L0, maxlevel, minlevel, rho1, rho2, mu1, mu2, 
                  f.sigma, gravity, ud0, f.wetting.theta_s,
                  xd0, df, ds, tend, tout, trestart, no_bubbles, DT);
}

event acceleration (i++) 
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= gravity;
}

event logfile (i++; t <= tend)
{
  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += f[]*sq(Delta)*cm[];

  if (i == 1) v0 = vreal;
  double verr = i == 0? 0: (vreal - v0)/v0 * 100;

  /**
  Calculate height of film at impact location (hf). */
  double hf = 0;
  foreach_boundary(bottom, reduction(max:hf)) {
    if (f[] > 1e-6 && f[] < 1-1e-6 && cs[] == 1) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double htemp = x + p.x*Delta - 0.5*ds;
      if (htemp > hf)
        hf = htemp;
    }
  }

  /** 
  Find contact line cell and calculate volume, thetacl, dmax, and ucl.
  Also calculate the total normal force acting on the sphere. (fn) */

  double betacl = 0, thetacl = 0, dmax = 0, dmax_abs = 0, fn = 0;

  foreach(reduction(max:betacl) reduction(max:thetacl) reduction(max:dmax) 
          reduction(max:dmax_abs) reduction(+:fn)) {

    if (extra[] || (on_interface(cs) && f[] > INT_TOL && f[] < cs[] - INT_TOL)) { // three-phase cell

      coord nf = interface_normal(point, f), ns = facet_normal(point, cs, fs), ts = {-ns.y, ns.x};
      double alphas = plane_alpha (cs[], ns), alphaf = plane_alpha (f[], nf);

      if (f.wetting.theta_s < 180) {

        normalize2(&ns); normalize2(&ts);
        nf = normal_contact (ns, nf, contact_angle[]);
        normalize_sum(&ns); normalize_sum(&nf);
        alphaf = immersed_alpha (f[], cs[], nf, 0, ns, alphas, f[]);

      }

      coord pint = {0,0};
      plane_area_center (nf, alphaf, &pint);
      if (extra[]) {
        int count = interface_intersect (nf, alphaf, ns, alphas, pint = &pint, 
                                         lhs = (coord){-0.6,-0.6}, rhs = (coord){0.6, 0.6});
        if (!count)
            continue;

        coord pintg = local_to_global_coord(point, pint);
        double betacl0 = atan2(pintg.y,pintg.x);

        double dmax0 = ds * betacl0; // total arc length s = 2 * (ds/2) * betacl
        if (dmax0 > dmax)
            dmax = dmax0;

        betacl0 *= 180./pi;

        if (betacl0 > betacl) {
            betacl = betacl0;
            thetacl = contact_angle[]*180./pi;
       	}
      }
    }

    // calculate the (normal) force acting on the sphere
    if (on_interface(cs)) {
        coord n = facet_normal(point, cs, fs), pint;
        foreach_dimension()
            n.x *= -1;
        double alpha = plane_alpha (cs[], n);
        double area = Delta*2*M_PI*y*plane_area_center(n, alpha, &pint);
        pint = local_to_global_coord(point, pint);
        double ps = interpolate_linear (point, p, pint.x, pint.y, 0);
        foreach_dimension()
          fn += ps*area*n.x;
    }

    // calculate the absolute maximum diameter
    if (f[] > 0 && f[] < 1) {
      coord nf = interface_normal(point, f), pint = {0,0};
      double alphaf = plane_alpha (f[], nf);
      plane_area_center (nf, alphaf, &pint);
      coord pintg = local_to_global_coord(point, pint);
      double dmax_abs0 = ds * atan2(pintg.y, pintg.x);
      if (dmax_abs0 > dmax_abs)
        dmax_abs = dmax_abs0;
    }
  }

  fprintf(stdout, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
    i, t, ds, df, f.wetting.theta_r, f.wetting.theta_a, v0, vreal, verr, hf, thetacl, betacl, dmax, dmax_abs, fn);
  fflush(stdout);
}

/**
Output dump files of the simulation. Here we dump the pressure in case we want
to make movies that involve it in the future. */
event dumps (t += tout)
{
  p.nodump = false;
  scalar * dumplist = {u,p,f,cs,ch,contact_angle,extra};
  char name[80];
  sprintf(name, "dumps/dump-%g", t);
  dump (file = name, list = dumplist);
}

event remove_bubbles (i+=50)
{
  if (no_bubbles) {
    double dc = 0;
    double dvol = 0;
    int count = 0;
    foreach(reduction(+:dc) reduction(+:dvol) reduction(+:count)) {
      if (f[] > 0 && f[] < cs[] && cs[] > 0 && cs[] < 1) {
        bool out = false;
        foreach_neighbor() {
          if (cs[] && !f[])
            out = true;
        }
        if (!out) {
          dc += cs[] - f[];
          dvol += dc*dv();
          f[] = cs[];
          count++;
        }
      }
    }
    if (dvol > 0)
      fprintf(stderr, "fixed %d bubbles. dc = %g dvol = %g\n", count, dc, dvol);
    boundary({f});
  
    remove_droplets(f, 5, bubbles=true);
  }
}

/**
For high We or on coarse grids, secondary droplets can form, either from physical
or numerical breakup. These droplets, if tiny, can cause problems, so we remove
them without impacting the overall dynamics. */
event remove_droplets (i += 50)
{
  remove_droplets(f, 5);
}

#if 0
scalar dis[];
event energy_budget (i++)
{
  /** initial energy */
  const double ke0 = 0.5*rho1*pi*cube(df)*sq(ud0)/6.;
  const double se0 = f.sigma*pi*sq(df);

  double ked = 0., keg = 0., sed = 0., sedw = 0., ket = 0., ped = 0., vdt = 0.;
  foreach(reduction(+:ket) reduction(+:sed) reduction(+:sedw)
          reduction(+:ked) reduction(+:ped) reduction(+:keg)
	  reduction(+:vdt)) {

    dis[] = 0;

    if (cs[] > 0.) {

      double dv_axi = pow(Delta, 2.)*2.*pi*y;

#if IBM
      double val = gc[];
#else
      double val = cs[] > 0.5? 1: 0;
#endif

      /** kinetic energy */
      double ke = 0.;
      foreach_dimension() {
        ke += sq(u.x[]);
      }
      ked += val*f[]*ke*rhov[]/(cs[])*dv_axi;   // droplet KE
      keg += val*(1-f[])*ke*rhov[]/(cs[])*dv_axi;   // droplet KE
      ket += val*ke*rhov[]/(cs[])*dv_axi;       // total KE

      /** surface energy */
      if (cs[] == 1 && f[] > 0 && f[] < 1) { // liquid-gas interface
        coord nf = interface_normal (point, f), p;
        double alphaf = plane_alpha(f[], nf);
        sed += pow(Delta, dimension - 1)*plane_area_center (nf, alphaf, &p)*2.*pi*y;
      }
      if (cs[] > 0 && cs[] < 1 && f[]) { // liquid-solid interface
        coord ns = facet_normal (point, cs, fs), p;
        double alphas = plane_alpha(cs[], ns);
        sedw += pow(Delta, dimension - 1)*plane_area_center (ns, alphas, &p)*2.*pi*y;
      }
      
      /** potential energy */
      ped += f[] * (rhov[]/cs[]) * gravity * (x - 0.5*ds) * dv_axi; // - 0.5*ds to make PE = 0 at the top of the sphere

      /** viscous dissipation */   
      double dvdr = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
      double dudx = (u.x[1] - u.x[-1])/(2.*Delta);
      double dvdx = (u.y[1] - u.y[-1])/(2.*Delta);
      double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);

      double dp = val*2.*mu(f[])*(sq(dvdr) + sq(u.y[]/y) + sq(dudx) + 0.5*sq(dvdx + dudy));
      
      dis[] = dp * dv_axi;
      vdt += dt * dp * dv_axi;
    }
  }

  ked *= 0.5;
  ket *= 0.5;
  keg *= 0.5;
  sed *= f.sigma;
  sedw *= -f.sigma * cos(f.wetting.theta_s*pi/180.);

  static FILE * fp = fopen("energy", "w");
  fprintf(fp, "%d %g %g %g %g %g %g %g %g %g %g\n",
               i, t, ke0, se0, ked, keg, ket, ped, sed, sedw, vdt);
  fflush(fp);
}
#endif

void print_film_profile(double t)
{
  char name[80];
  sprintf(name, "data/%d-film-%g", pid(), t);
  FILE * fp = fopen(name, "w");

  vector nf[];
  scalar alphaf[];
  reconstruction(ch, nf, alphaf);

  foreach() {
    if (cs[] == 1 && f[] > 1e-4 && f[] < 1 - 1e-4) {
      coord nf0 = {nf.x[], nf.y[]}, p;
      plane_area_center (nf0, alphaf[], &p);
      p = local_to_global_coord(point, p);

      double hf = magnitude_coord(p) - 0.5*ds;
      double angle = atan2(p.y, p.x);

      fprintf(fp, "%g %g %g %g\n", p.x, p.y, hf, angle);
    }
  }
  fclose(fp);
}

event film_profile (t = {0.1896, 0.473, 0.757, 1.0365, 1.1804})
{
  print_film_profile(t);
}

event movie (t += tout) 
{
  /** calculate vorticity only for fluid cells */
  scalar omega[];
  vorticity2(u, omega);

  char name[80];

  double fov0 = 2*2*atan(3*df/2.);
  double tx0 = -0.5*ds/L0;

  /** calculate min and max */
  stats pn = statsf(p);
  stats on = statsf(omega);

  double pbound = pn.max;
  if (fabs(pn.max) < fabs(pn.min))
    pbound = fabs(pn.min);

  double obound = on.max;
  if (fabs(on.max) < fabs(on.min))
    obound = fabs(on.min);

  /** identify furthest contact line cell */
  double xi = 0;
  foreach(reduction(max:xi)) {
    if (extra[]) {
        double xi0 = atan2(y,x);
        if (xi0 > xi)
          xi = xi0;
    }
  }

  coord pcl = {0,0};
  foreach(serial) {
    if (extra[] && atan2(y,x) == xi) {
      coord n = facet_normal (point, cs, fs), pint;
      double alpha = plane_alpha (cs[], n);
      plane_area_center (n, alpha, &pint);
      pcl = local_to_global_coord(point, pint);
    }
  }

  /** vorticity */
  view(fov = fov0, tx = tx0, width = 1000, height = 1000);
  clear();
  draw_vof ("ch", lw = 2);
  draw_vof("cs", "fs");
  squares("ch", min = -1, max = 1);
  cells();
  mirror (n = {0,-1}) {
    draw_vof("ch", lw = 2);
    draw_vof("cs", "fs");
  squares("omega", min = -obound/20, max = obound/20, cbar = true, map = cool_warm);
  }
  sprintf(name, "imgs/vort-%.5f.png", t);
  save(name);

  /** pressure */
  view(fov = fov0, tx = tx0, width = 1000, height = 1000); 
  clear();
  draw_vof ("ch", lw = 2);
  draw_vof("cs", "fs");
  squares("ch", min = -1, max = 1);
  mirror (n = {0,-1}) {
    draw_vof("ch", lw = 2);
    draw_vof("cs", "fs");
    squares("p", min = -pbound/5, max = pbound/5, cbar = true, map = cool_warm);
  }
  sprintf(name, "imgs/pres-%.5f.png", t);
  save(name);

  /** contact line */
  view(fov = 0.1, tx = -pcl.x/L0, ty = -pcl.y/L0, width = 1000, height = 1000);
  //view(fov = 30, tx = -pcl.x/L0, ty = -pcl.y/L0, width = 1000, height = 1000, far=1000, near=1, tz=-0.01);
  clear();
  draw_vof ("ch", lw = 2);
  //squares ("ch");
  draw_vof("cs", "fs", lw = 2);
  cells();
  vectors("uf", scale=5e-4, lw = 1.5);
  labels("extra");
  sprintf(name, "imgs/cl-%.5f.png", t);
  save(name);

  print_film_profile(t);
}


#if 0

scalar fpos[];
scalar npos[];

/** output data relating to the film and rim */
event rim_data (i+=5)
{
  static FILE * fp = fopen ("rim", "w");

  vector nf[];
  scalar alphaf[];
  reconstruction(ch, nf, alphaf);

  /** film statistics */
  double favg = 0;
  int count = 0;
  foreach(reduction(+:favg) reduction(+:count)) {
    if (cs[] == 1 && ch[] > 0 && ch[] < 1) {
      coord nf0 = {nf.x[], nf.y[]}, p;
      plane_area_center (nf0, alphaf[], &p);
      p = local_to_global_coord(point, p);

      double r = sqrt(sq(p.x) + sq(p.y)) - 0.5*ds;
      fpos[] = r;
      favg += fpos[];
      count++;

      coord nr = {p.x, p.y};
      coord tr = {-p.y, p.x};

      double dpt = fabs(dot_product_norm(nf0, tr));
      double dpn = fabs(dot_product_norm(nf0, nr));

      if (dpt > dpn)
        npos[] = nodata;
      else
        npos[] = dot_product_norm(nf0, nr);
    }
    else {
        fpos[] = nodata;
        npos[] = nodata;
    }
  }

  favg /= (double)(count);
  stats filmstat = statsf(fpos);

  /** rim statistics */
  const double hmin = L0 / (1 << maxlevel);
  double hrim = hmin*3, hrim2 = hmin*3; // calculate height of rim
  double anglerim = 0, anglerim2 = 0; // angle of max rim height
  double anglemax = 0; // angle of maximum film height
  foreach(reduction(max:anglemax) reduction(max:anglerim) reduction(max:anglerim2)) {
    if (fpos[] != nodata && y > 0.01*ds) { // y > 0.01*ds to avoid droplet on sym. boundary
      bool max = true;
      double fpos0 = fpos[];
      foreach_neighbor(1) {
        if ((fpos[] != nodata && fpos[] > fpos0) || (on_interface(cs))) {
          max = false;
          break;
        }
           
      }
      if (max) {
        coord nf0 = {nf.x[], nf.y[]}, p;
        plane_area_center (nf0, alphaf[], &p);
        p = local_to_global_coord(point, p);
        double anglerim0 = atan2(p.y,p.x);
        if (anglerim0 > anglerim) {
            anglerim = anglerim0;
        }
      }
    }
    if (npos[] != nodata && y > 0.01*ds) { // y > 0.01*ds to avoid droplet on sym. boundary
      bool max = true;
      double npos0 = npos[];
      foreach_neighbor(1) {
        if ((npos[] != nodata && npos[] > npos0) || (on_interface(cs))) {
          max = false;
          break;
        }
           
      }
      if (max) {
        coord nf0 = {nf.x[], nf.y[]}, p;
        plane_area_center (nf0, alphaf[], &p);
        p = local_to_global_coord(point, p);
        double anglerim0 = atan2(p.y,p.x);
        if (anglerim0 > anglerim2) {
            anglerim2 = anglerim0;
        }
      }
    }
    if (fpos[] == filmstat.max) {
      coord nf0 = {nf.x[], nf.y[]}, p;
      plane_area_center (nf0, alphaf[], &p);
      p = local_to_global_coord(point, p);
      anglemax = atan2(p.y,p.x);
    }
  }

  foreach(reduction(max:hrim) reduction(max:hrim2)) {
    if (fpos[] != nodata && y > 0.01*ds) {
      coord nf0 = {nf.x[], nf.y[]}, p;
      plane_area_center (nf0, alphaf[], &p);
      p = local_to_global_coord(point, p);
      double anglerim0 = atan2(p.y,p.x);

      if (anglerim0 == anglerim) {
        hrim = fpos[];
      }
    }
    if (npos[] != nodata && y > 0.01*ds) {
      coord nf0 = {nf.x[], nf.y[]}, p;
      plane_area_center (nf0, alphaf[], &p);
      p = local_to_global_coord(point, p);
      double anglerim0 = atan2(p.y,p.x);

      if (anglerim0 == anglerim2) {
        hrim2 = sqrt(sq(p.x) + sq(p.y)) - 0.5*ds;
      }
    }

  }
  
  fprintf(fp, "%d %g %g %g %g %g %g %g %g %g %g\n", 
    i, t, filmstat.min, filmstat.max, filmstat.stddev, favg, anglemax, hrim, anglerim, hrim2, anglerim2);
  fflush(fp);
}
#endif

/**
We use an intermediate field f1 to trick the adapt function into always
keeping the liquid interface at the maximum resolution. Note we do not do
this for the solid interface.

Any cell with fluid volume fraction (cs) < 1e-3 is removed and set to 0. */

event adapt (i++)
{
  scalar f1[];
  foreach() 
    f1[] = f[];

  adapt_wavelet ({cs, f1, u}, (double[]){1e-3, 1e-3, uemax, uemax},
    maxlevel = maxlevel, minlevel = minlevel);

  solid (cs, fs, (sq(x) + sq(y) - sq(ds*0.5)));

  int smalls = fractions_cleanup(cs, fs, smin = 1e-3);
  if (smalls > 0)
    fprintf(stderr, "changed %d small solid cells\n", smalls);
}

