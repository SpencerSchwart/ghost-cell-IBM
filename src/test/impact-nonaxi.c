#define CA 1
#include "grid/octree.h"
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../contact-ibm.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#if TRACE
#include "profiling.h"
#endif

#define MAXLEVEL 8
#define MINLEVEL 5

double xc = 0;
double yc = 0;
double zc = 0;

double hi0;             // inital droplet position
double ud0 = 1;         // initial drop velocity
double df = 1;          // drop diameter
double ds = 2.727;       // sphere diameter
double gravity = 0.011499;

double tmovie = 0.0125;
double tdump = 0.1;
double t_end = 12;
double theta0 = 30;

double alpha0 = 36;     // angle of attack
double alpha0_tol = 0.5; // +/- degree window for calculating h_i

u.n[immersed] = dirichlet(0);
u.t[immersed] = dirichlet(0);
u.r[immersed] = dirichlet(0);

p  [right] = dirichlet(0);
pf [right] = dirichlet(0);
u.n[right] = neumann(0);

p  [left] = dirichlet(0);
pf [left] = dirichlet(0);
u.n[left] = neumann(0);

u.n[back] = dirichlet(0);
f[back] = neumann(0);

int main() {
  size(10);

  origin (-10./4., -10./2., 0);

  init_grid (1 << MINLEVEL);
  
  rho1 = 1;
  rho2 = 0.001532;
  mu1 = 1.02E-3;
  mu2 = 7.55E-6;
  f.sigma = 0.00647;

  TOLERANCE = 1e-4;
  
  const scalar c[] = theta0*pi/180.;
  contact_angle = c;
  run();
}

event init (t = 0)
{
  hi0 = (0.5*ds + 0.5*df + 0.05); // initial distance of droplet
  double alpha = alpha0*pi/180.;

  coord ns = {cos(alpha), sin(alpha), 0};
  
  coord x0; // inital droplet position
  foreach_dimension()
    x0.x = hi0*ns.x;

  scalar ibm1[], f1[];
  face vector ibmf1[];

  int count = 0;
  astats st;
  do {
    solid (ibm1, ibmf1, (sq(x - xc) + sq(y - yc) + sq(z - zc) - sq(ds*0.5)));
    fraction (f1, - (sq(x - x0.x) + sq(y - x0.y) + sq(z - x0.z) - sq(df*0.5)));
    st = adapt_wavelet ({ibm1, f1}, (double[]){1e-5, 1e-5}, MAXLEVEL, MINLEVEL);
    count++;
  } while ((st.nc || st.nf) && count < 40);

  solid (ibm, ibmf, (sq(x - xc) + sq(y - yc) + sq(z - zc) - sq(ds*0.5)));
  fractions_cleanup(ibm, ibmf, 1e-4);

  fraction(f,  - (sq(x - x0.x) + sq(y - x0.y) + sq(z - x0.z) - sq(df*0.5)));
  fraction(ch, - (sq(x - x0.x) + sq(y - x0.y) + sq(z - x0.z) - sq(df*0.5)));
  fraction(cr, - (sq(x - x0.x) + sq(y - x0.y) + sq(z - x0.z) - sq(df*0.5)));

  foreach() {
    u.x[] = -ud0;
    //u.x[] = -(ud0*cos(alpha))*f[];
    //u.y[] = -(ud0*sin(alpha))*f[];
    //u.z[] = 0;
  }
}

event acceleration (i++) 
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= gravity;
}

double v0 = 0;
event logfile (i++; t <= t_end)
{
  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += cr[]*cube(Delta);

  if (i == 1) v0 = vreal;

  double alpha = alpha0*pi/180.;
  coord ns = {cos(alpha), sin(alpha), 0};

  double hf = -HUGE;
  foreach_boundary(back, reduction(max:hf)) {
    if (ch[] > 1e-3 && ch[] < 1 && ibm[] && cr[] > 1e-3) {
        coord n = interface_normal (point, ch), p;
        double alpha = plane_alpha (ch[], n);
        plane_area_center (n, alpha, &p);
        coord cc = {x,y,z}, pp;
        foreach_dimension()
            pp.x = cc.x + p.x*Delta;

        double h = pp.x*ns.x + pp.y*ns.y + pp.z*ns.z;
        double apparent_angle = atan(pp.y/pp.x)*180./pi;
        if (h > hf && fabs(apparent_angle) > alpha0 - alpha0_tol)
            hf = h;
    }
  }

  double verr = i == 0? 0: (vreal - v0)/v0 * 100;

  fprintf(stderr, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verr, hf/hi0);
  fprintf(stdout, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, verr, hf/hi0);
}

event dump (t += tdump)
{
  scalar * dumplist = {u,p,f,ibm,cr,ch};
  char name[80];
  sprintf(name, "dump-%g", t);
  dump (file = name, list = dumplist);
}

event movie (t += tmovie) 
{
    char name[80];
    sprintf(name, "%g-front-%g.png", theta0, t);
    view(camera="front", fov = 11, tx = -0.075, width = 1000, height = 1000, bg={1,1,1});
    draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
    draw_vof("ibm", "ibmf", fc={0.796875, 0.796875, 0.796875});
    save(name);

    sprintf(name, "%g-back-%g.png", theta0, t);
    view(camera="back", fov = 11, tx = -0.075, width = 1000, height = 1000, bg={1,1,1});
    draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
    draw_vof("ibm", "ibmf", fc={0.796875, 0.796875, 0.796875});
    save(name);

    float q[4] = {0, 0, 0, 0};
    gl_axis_to_quat((float[]){0,1,0},pi/4., q);
    sprintf(name, "%g-iso-%g.png", theta0, t);
    view(quat=q, fov = 11, tx = -0.075, width = 1000, height = 1000, bg = {1,1,1});
    draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
    draw_vof("ibm", "ibmf", fc={0.796875, 0.796875, 0.796875});
    mirror (n = {0,0,-1}) {
      draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
      draw_vof("ibm", "ibmf", fc={0.796875, 0.796875, 0.796875});
    }
    save(name);

    sprintf(name, "%g-top-%g.png", theta0, t);
    view(camera="top", fov = 11, tx = -0.075, width = 1000, height = 1000, bg = {1,1,1});
    draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
    draw_vof("ibm", "ibmf", fc={0.796875, 0.796875, 0.796875});
    mirror (n = {0,0,-1}) {
      draw_vof ("ch", lw = 2, fc={1, 0.62891, 0.62891});
      draw_vof("ibm", "ibmf", fc={0.796875, 0.796875, 0.796875});
    }
    save(name);
}

event adapt (i++)
{
  scalar f1[], ibm1[];
  foreach() {
    f1[] = ch[];
    ibm1[] = ibm[];
  }
  adapt_wavelet ({ibm1, f1, u}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-3}, 
    maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

