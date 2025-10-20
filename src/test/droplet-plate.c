#define CA 1
#define MOVING 1
#define MOVIE 0
#define SHAPE 1

#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
//#include "../contact-ibm.h"
#include "../contact-ibm3D.h"
#include "view.h"

#define LEVEL 8
#define MIN_LEVEL 4
#define L0 4.
#define D  0.5             // diameter of drop
#define LP 5               // length of plate
#define HP 0.1              // height of plate

const double theta0 = 75;
//const double t_end = 0.50001;
const double t_end = 0.70001;

//#define upy(t) (9.81*(t))
#define upy(t) (1)

//coord fi = {1., 2./10. + D/3.}; // initial drop position
coord fi = {1., 2./10.};

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

u.t[left] = neumann(0);
u.n[left] = neumann(0);

u.n[top] = neumann(0);
u.t[top] = neumann(0);

u.n[bottom] = neumann(0);
u.t[bottom] = neumann(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);

void plate_geometry (scalar ibm, face vector ibmf, coord xp = {0,0})
{
    vertex scalar phi[];
    foreach_vertex() {
        double a = min(x + 0.5*LP, -x + 0.5*LP);
        double b = min((y - xp.y) - 2./10., HP + 2./10. - (y - xp.y));

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
        if (y > 2./10. + HP/2.)
            phi[] = -sq(x + fi.x) - sq(y - fi.y) + sq(D/2.);
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

    origin (-L0/2, 0);

    //rho1 = 1.225; rho2 = 0.1694;
    //mu1 = mu2 = 0.00313;
    rho1 = rho2 = 1;
    mu1 = mu2 = 0.1;

    f.sigma = 1.;
    
    TOLERANCE = 1e-6;

#if !MOVING 
    a[] = {0, -9.81};
#endif

    CFL = 0.5;

    const scalar c[] = theta0*pi/180;
    contact_angle = c;
    run();
}

double v0 = 0;
event init (t = 0)
{
    mask (x > 0? right: none);

    plate_geometry (ibm, ibmf);
    event ("update_metric");

    droplet_geometry (f);
    
    v0 = real_volume(f);
}

#if MOVING
void move_solid_y (scalar ibm, face vector ibmf)
{
    //coord xp = {0, (sq(t + dt/2.))/2. * 9.81};
    coord xp = {0, (t + dt/2.) * upy(t)};
    //fprintf (stderr, "moving plate @t=%0.15g dt=%0.15g t_s=%0.15g w/u=%0.15g xp.y=%0.15g\n",
    //    t, dt, t_start, upy(t), xp.y);
    plate_geometry (ibm, ibmf, xp);
    event ("update_metric");
}

void move_solid_x (scalar ibm, face vector ibmf)
{
}
#endif

event logfile (i++; t <= t_end)
{
  
  double vreal = real_volume(f);

  coord xp = {0, (sq(t + dt/2.))/2. * 9.81};

  fprintf(stderr, "%d %g %g %g %g %g %g %g %g\n", 
    i, t, dt, theta0, v0, vreal, (vreal - v0)/(v0+SEPS)*100, upy(t), xp.y);
  fprintf(stdout, "%d %g %g %g %g %g %g %g %g\n", 
    i, t, dt, theta0, v0, vreal, (vreal - v0)/(v0+SEPS)*100, upy(t), xp.y);
}

#if MOVIE
event movie (i += 10) 
{
    #if 1
    char name[80];
    sprintf (name, "movie-%d.png", i);
    FILE * fp1 = fopen (name, "w");
    view (fov = 20, camera = "front", tx = 0., ty = -0.5);
    clear();
    draw_vof ("f", lw = 2);
    draw_vof ("ibm", "ibmf", filled=-1);
    squares("f", min = 0, max = 1);
    save(fp = fp1);
    fclose (fp1);
    #endif
}
#endif

#if SHAPE
event interface (i += 100)
//event interface (t += 0.1)
{
    char name[80];
    sprintf (name, "%d-shape-%d", MOVING, i);
    FILE * fp = fopen (name, "w");
    output_facets (f, fp);
    fclose (fp);

    sprintf (name, "%d-ibmshape-%d", MOVING, i);
    FILE * fp1 = fopen (name, "w");
    output_facets (ibm, fp1);
    fclose (fp1);
}

event interface_end (t = end)
{
    char name[80];
    sprintf (name, "%d-shape-end", MOVING);
    FILE * fp = fopen (name, "w");
    output_facets (f, fp);
    fclose (fp);

    sprintf (name, "%d-ibmshape-end", MOVING);
    FILE * fp1 = fopen (name, "w");
    output_facets (ibm, fp1);
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
