#define CA 1

#include "../ibm-gcm.h"
#include "../my-axi.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
#include "../ibm-gcm-vof-test.h"
#include "../my-two-phase.h"
#include "../my-tension.h"
//#include "../my-conserving.h"
#include "../contact-ibm.h"

#define D0 1.

u.t[immersed] = neumann(0);
u.n[immersed] = dirichlet(0);

u.n[top] = dirichlet(0);

p[right] = dirichlet(0);
p[right] = dirichlet(0);
u.n[right] = neumann(0);

int max_level = 9;
int min_level = 5;
double L = 4.;
double t_end = 100;
double theta0 = 90.01;
double amp = 0.e-2;
double dia = 5.e-3;
double gpert = 4.9;

double femax = 0.01;
double semax = 0.01;
double uemax = 0.001;

double rhol  = 1000;
double rhog  = 1.2;
double mul   = 1.e-3;
double mug   = 1.8e-5;
double sigma = 0.07;

double omegac = 0.5367;
double time_pert = 1.0434;

double Oh = 0.005;

int main()
{
  size (L);
  origin (-0.2626,0); // small cell @lvl8
  // origin (-0.26, 0); // no small cell @lvl8

  init_grid (1 << max_level);
  
  TOLERANCE = 1e-5;
  rho1 = rhol/rhol;
  rho2 = rhog/rhol;
  //mu1 = mul/sqrt(rhol*dia*sigma);
  mu1 = Oh*sqrt(rhol*dia*sigma);
  mu2 = mug/sqrt(rhol*dia*sigma);
  f.sigma = sigma/sigma;

  const scalar c[] = theta0*pi/180.;
  contact_angle = c;

  run();
}

double v0 = 0;

event init (t = 0)
{

  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  
  vertex scalar phi[];
  foreach_vertex()
    phi[] = x;
  boundary ({phi});
  fractions (phi, ibm, ibmf);
  fraction (f, - (sq(x) + sq(y) - sq(D0/2.)));
  fraction (ch, - (sq(x) + sq(y) - sq(D0/2.)));

  v0 = real_volume(f);
}

event logfile (i++; t <= t_end)
{
  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += cr[]*sq(Delta)*cm[];

  scalar pos[];
  position (ch, pos, {1,0});
  double hmax = statsf(pos).max;

  double perror = i == 0? 0: (vreal - v0)/v0 * 100;

  fprintf(stderr, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, perror, hmax);
  fprintf(stdout, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, perror, hmax);
}


event acceleration (i++)
{
  foreach_face(x)
    a.x[] += (t > time_pert? 0:
               -gpert*rhol*sq(dia)/sigma);
}


event adapt (i++)
{
  scalar ch1[], ibm1[];
  foreach()
    ch1[] = ch[], ibm1[] = ibm[];
    
  adapt_wavelet ({ch1,ibm1,u}, (double[]){femax,semax,uemax,uemax}, 
                  minlevel = min_level, maxlevel = max_level);
}

