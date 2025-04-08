#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
#include "../contact-ibm.h"

#define LEVEL 7
#define L0 2.
#define D  0.25             // diameter of drop
#define LP 1.               // length of plate
#define HP 0.1              // height of plate

const double theta0 = 90.;
const double t_end = 10.;
const double t_start = 0.25;

coord fi = {LP, L0/10. + HP + D/3.}; // initial drop position
coord vp = {0., 0.01};               // velocity of plate

// Boundar Conditions
u_x_ibm_dirichlet (0)
u_y_ibm_dirichlet (t >= t_start? vp.y: 0)

p[top] = dirichlet(0);
pf[top] = dirichlet(0);

u.n[bottom] = dirichlet(0);

void plate_geometry (scalar ibm, face vector ibmf, coord xp = {0})
{
    vertex scalar phi[];
    foreach_vertex() {
        double a = min(x - LP/2., 1.5*LP - x);
        double b = min((y - xp.y) - L0/10., HP + L0/10. - (y - xp.y));

        phi[] = max(-a,-b);
    }
    boundary ({phi});
    fractions (phi, ibm, ibmf);
    boundary({ibm, ibmf});
}

void droplet_geometry (scalar f)
{
    vertex scalar phi[];
    foreach_vertex() {
        if (y > L0/10. + HP/2.)
            phi[] = -sq(x - fi.x) - sq(y - fi.y) + sq(D/2.);
        else
            phi[] = -1;
    }
    fractions (phi, f);
    boundary({f});
}

int main()
{
    size (L0);
    init_grid (1 << LEVEL);

    mu1 = mu2 = 0.1;
    f.sigma = 1.;

    const scalar c[] = theta0*pi/180;
    contact_angle = c;
    run();
}

double v0 = 0;
event init (t = 0)
{
    plate_geometry (ibm, ibmf);
    droplet_geometry (f);
    
    v0 = real_volume(f);
}

#if 0
event tracer_advection (i++)
{
    coord xp = {0, t >= t_start? (t - t_start)*vp.y: 0};
    plate_geometry (ibm, ibmf, xp);
    event ("update_metric");
}
#endif

event logfile (i++; t <= t_end)
{
  double vreal = real_volume(f), vreal2 = 0;
  foreach()
    vreal2 += cr[]*sq(Delta);

  fprintf(stderr, "%d %g %g %g %g %g\n", i, t, theta0, v0, vreal, vreal2);  
}
