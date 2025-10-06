#define CA 1
#define RAIN 1

#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../contact-ibm.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "view.h"

#define LEVEL 6
#define MIN_LEVEL 4
#define L0 (10)
#define LP (L0)               // length of plate
#define HP 1            // height of plate

const double theta0 = 30;
const double t_end = 20;

// Boundary Conditions
u.n[immersed] = dirichlet(0);
u.t[immersed] = dirichlet(0);

p[top]      = dirichlet(0);
pf[top]     = dirichlet(0);

u.n[top]    = neumann(0);
u.t[top]    = neumann(0);

u.n[bottom] = neumann(0);
u.t[bottom] = neumann(0);

p[bottom]   = dirichlet(0);
pf[bottom]  = dirichlet(0);

u.n[left]   = dirichlet(0);
u.n[right]  = dirichlet(0);

void plate_geometry (scalar ibm, face vector ibmf, coord xp = {0,0})
{
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = y - 2.5;
    }
    boundary ({phi});
    fractions (phi, ibm, ibmf);
    boundary({ibm, ibmf});
}

void film_geometry (scalar f)
{
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = 5 - y < 0? -1: 1;
    }
    fractions (phi, f);
    boundary ({f});
}


int main()
{
    size (L0);
    init_grid (1 << LEVEL);

    origin (-L0/2, 0);

    rho2 = 1.225;
    rho1 = 0.1694;
    mu1 = mu2 = 0.00313;
    f.sigma = 1;

    // rain parameters
    f.urain = (coord){0, -1};
    f.rvf = 0.5;
    f.stddev = 0.5;
    f.mean = 0;
    f.normal = (coord){0,1};
    f.intercept = 0;

    TOLERANCE = 1e-6;

    const scalar c[] = theta0*pi/180;
    contact_angle = c;
    run();
}

double v0 = 0;
event init (t = 0)
{
    plate_geometry (ibm, ibmf);
    event ("update_metric");

    film_geometry (f);

}

event logfile (i++; t <= t_end)
{
 
  double vreal = 0;
  foreach()
    vreal += f[]*sq(Delta);

  if (i == 1) v0 = vreal;

  scalar pos[];
  position (f, pos, {0,1});
  double hmax = statsf(pos).max;

  fprintf(stderr, "%d %g %g %g %g %g %g %g\n", 
    i, t, dt, theta0, v0, vreal, (vreal - v0)/(v0+SEPS)*100, hmax);
  fprintf(stdout, "%d %g %g %g %g %g %g %g\n", 
    i, t, dt, theta0, v0, vreal, (vreal - v0)/(v0+SEPS)*100, hmax);
}

#if 0
//event movie (t += 0.11) // messes up timestep which causes errors
event movie (i += 50) 
{
    //double tyy = -0.5 - ((sq(t + dt/2.))/2. * 9.81) / L0;
    char name[80];
    sprintf (name, "movie-%d.png", i);
    FILE * fp1 = fopen (name, "w");
    view (fov = 20, camera = "front", tx = 0., ty = -0.5); // fov 20 for overall domain
    clear();
    draw_vof ("f", lw = 2);
    draw_vof ("ibm", "ibmf", filled=-1);
    squares("f", min = 0, max = 1);
    save(fp = fp1);
    fclose (fp1);
}
#endif

#if 1
event adapt (i++)
{
    scalar f1[], ibm1[];

    foreach() {
        f1[] = ch[];
        ibm1[] = ibm[];
    }
    adapt_wavelet ({ibm1, f1}, (double[]){1e-5,1e-5}, 
                    maxlevel = LEVEL, minlevel = MIN_LEVEL);
}
#endif
