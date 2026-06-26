#include "ibm/src/ibm-gcm.h"

scalar p[], pf[];
vector u[];
face vector uf[];

attribute {
  double sigma;
}

double mu1;

#include "run.h"
#include "timestep.h"

event defaults (i = 0)
{

  CFL = 0.2;

  uf.x.refine = refine_face;
  foreach_dimension() {
    uf.x.prolongation = refine_ibm_face_x;
  }
  for (scalar s in {u}) {
    s.restriction = restriction_ibm_linear;
    s.refine = s.prolongation = refine_ibm_linear;
    s.depends = list_add (s.depends, cs);
  }
}

double dtmax;

double t0 = 1. [0, 1];

event init (i = 0)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  
  dtmax = DT;
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (timestep (uf, dtmax));
}

event vof (i++,last);

scalar f[], * interfaces = {f};

#include "ibm/src/contact-ibm.h"
#include "ibm/src/my-vof-ca.h"

double radius, distance;
double theta0 = 150;

int main()
{
  origin (-0.5, -0.5);
  init_grid (1 << 7);

  f.wetting.theta_s = theta0;

  radius = 0.1409;
  distance = (sqrt(2*(1-cos(pi/6.)))*radius);

  run();
}

event init (i = 0)
{
  solid (cs, fs, (sq(x) + sq(y) - sq(radius)));
  fraction (f, - (sq(x) + sq(y-distance) - sq(radius)));
  fraction (ch, - (sq(x) + sq(y-distance) - sq(radius)));

    cs.dirty = true;
    foreach() {
        if (cs[] > GCV) // fluid cell
            gc[] = 1;
        else // ghost or solid cell
            gc[] = 0;
    }

    foreach_face() {
       if (cs[] > 0 || cs[-1] > 0 || is_ghost_cell (point, cs))
            gcf.x[] = 1;
       else
            gcf.x[] = 0;
    }
    boundary({gcf, gc});

  foreach() {
    u.x[] =  2.*pi*y*gc[];
    u.y[] = -2.*pi*x*gc[];
  }
}

event logfile (i++; t <= 1) {
  foreach_face()
    if (!fs.x[])
      uf.x[] = 0;
}

