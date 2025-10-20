#define CA 1
#define MOVING 0
#define MOVIE 0
#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof.h"
#include "../my-two-phase.h"
//#include "../contact-ibm.h"
#include "../contact-ibm3D.h"
#include "view.h"

#define LEVEL 9
#define MIN_LEVEL 4
#define L0 (8*M_PI)
#define LP (L0)               // length of plate
#define HP 1              // height of plate
#define RHOR (5.)

const double theta0 = 30;
const double t_end = 1.4;

#define upy(t) (9.81*(t))

// Boundary Conditions
u_x_ibm_dirichlet (0)

#if MOVING
u_y_ibm_dirichlet (upy(tt + dt))
#else
u_y_ibm_dirichlet (0)
#endif

p[top] = dirichlet(0);
pf[top] = dirichlet(0);

p[bottom] = dirichlet(0);
pf[bottom] = dirichlet(0);

u.n[left] = dirichlet(0);
uf.n[left] = dirichlet(0);

u.n[right] = dirichlet(0);
uf.n[right] = dirichlet(0);

u.n[top] = neumann(0);
u.t[top] = neumann(0);

u.n[bottom] = neumann(0);
u.t[bottom] = neumann(0);

void plate_geometry (scalar ibm, face vector ibmf, coord xp = {0,0})
{
    vertex scalar phi[];
    foreach_vertex() {
        double a = min(x + 0.5*LP, -x + 0.5*LP);
        double b = min((y - xp.y) - 0.6, HP - (y - xp.y));

        phi[] = max(-a,-b);
    }
    boundary ({phi});
    fractions (phi, ibm, ibmf);
    boundary({ibm, ibmf});
}

void film_geometry (scalar f)
{
    vertex scalar phi[];
    foreach_vertex() {
        if (y > 0.8)
            phi[] = (-cos(2*x)/2. + L0/2. - y);
        else
            phi[] = -1;
    }
    fractions (phi, f);
    boundary ({f});
}

int main()
{
    size (L0);
    init_grid (1 << LEVEL);

    origin (-L0/2, 0);

#if !MOVING
    a[] = {0,-9.81}; 
#endif

    rho2 = 1.225;
    rho1 = 0.1694;
    mu1 = mu2 = 0.00313;

    TOLERANCE = 1e-6;
    DT = 5e-4;
    CFL = 0.5;

    const scalar c[] = theta0*pi/180;
    contact_angle = c;
    run();
}

double v0 = 0;
event init (t = 0)
{
    mask (x > pi/2.? right: x < -pi/2.? left: none);
    plate_geometry (ibm, ibmf);
    event ("update_metric");

    film_geometry (f);
    
    v0 = real_volume(f);
}

#if MOVING
void move_solid_y (scalar ibm, face vector ibmf)
{
    coord xp = {0, (sq(t + dt/2.))/2. * 9.81};
    fprintf (stderr, "moving plate @t=%0.15g dt=%0.15g w/u=%0.15g xp.y=%0.15g\n",
        t, dt, upy(t), xp.y);
    plate_geometry (ibm, ibmf, xp);
    event ("update_metric");
}

void move_solid_x (scalar ibm, face vector ibmf)
{
}
#endif

scalar rhot[];

event logfile (i++; t <= t_end)
{
 
  foreach()
    rhot[] = rho[];

  double vreal = real_volume(f);

  fprintf(stderr, "%d %g %g %g %g %g %g %g\n", 
    i, t, dt, theta0, v0, vreal, (vreal - v0)/(v0+SEPS)*100, upy(t));
  fprintf(stdout, "%d %g %g %g %g %g %g %g\n", 
    i, t, dt, theta0, v0, vreal, (vreal - v0)/(v0+SEPS)*100, upy(t));
}

#if MOVIE
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
        f1[] = vertex_average(point, f);
        ibm1[] = ibm[];
    }
    adapt_wavelet ({ibm1, f1}, (double[]){1e-5,1e-15}, 
                    maxlevel = LEVEL, minlevel = MIN_LEVEL);
}
#endif
