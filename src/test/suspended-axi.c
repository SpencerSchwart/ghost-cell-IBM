scalar f[], * interfaces = {f};

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-generic.h"
#include "navier-stokes/conserving.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"

#define D0 1.

vector h[];
const double theta0 = 90;
h.t[left] = contact_angle ( theta0*pi/180.);

u.n[bottom] = dirichlet(0);

u.n[top] = dirichlet(0);

p[right] = dirichlet(0);
p[right] = dirichlet(0);
u.n[right] = neumann(0);


int max_level = 7;
int min_level = 5;
double L = 2.;
double t_end = 100;
double amp = 0.e-2;
double dia = 5.e-3;
double gpert = 4.9;

double femax = 0.01;
double uemax = 0.001;

double rhol  = 1000;
double rhog  = 1.2;
double mul   = 1.e-3;
double mug   = 1.8e-5;
double sigma = 0.07;

double omegac = 0.5367;
double time_pert = 1.0434;


int main()
{
  size (L);
  init_grid (1 << max_level);
  
  TOLERANCE = 1e-5;
  rho1 = rhol/rhol;
  rho2 = rhog/rhol;
  mu1 = mul/sqrt(rhol*dia*sigma);
  mu2 = mug/sqrt(rhol*dia*sigma);

  f.sigma = sigma/sigma;

  f.height = h;

  run();
}

double v0 = 0;
event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(D0/2.)));
  v0 = statsf(f).sum;
}

#if 1
event acceleration (i++)
{
  foreach_face(x)
    a.x[] += (t > time_pert? 0:
               -gpert*rhol*sq(dia)/sigma);
}
#endif

event logfile (i++; t <= t_end)
{
  double vreal = 0;
  foreach(reduction(+:vreal))
    vreal += f[]*sq(Delta)*cm[];

  scalar pos[];
  position (f, pos, {1,0});
  double hmax = statsf(pos).max;

  double perror = i == 0? 0: (vreal - v0)/v0 * 100;

  fprintf(stderr, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, perror, hmax);
  fprintf(stdout, "%d %g %g %g %g %g %g\n", i, t, theta0, v0, vreal, perror, hmax);
}


event adapt (i++)
{
  scalar f1[];
  foreach()
    f1[] = f[];
    
  adapt_wavelet ({f1,u}, (double[]){femax,uemax,uemax}, 
                  minlevel = min_level, maxlevel = max_level);
}

