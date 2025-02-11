#include "../ibm-gcm.h"
#include "../my-centered.h"
#include "../ibm-gcm-events.h"
//#include "navier-stokes/double-projection.h"
#include "view.h"

#define MAX_THICKNESS 0.15
#define CHORD_LENGTH 1
#define L0 20.
#define Re (10000.)
#define LEVEL 15
#define MIN_LEVEL 7

const double U0 =  1.0; // inlet velocity
const double rr = 1.1019*sq(MAX_THICKNESS); // Radius of leading edge
const double trailingOffset = 0.98;
const double theta_p = 0*pi/180; // initial A.O.A = 0 degrees
const double pitchStartTime = 0.2;
const double rotationRate = 0.6;
const double t_end = 2.;

coord ci = {L0/4., L0/2.}; // initial coordinates of airfoil
coord cr = {0.25*(CHORD_LENGTH), 0.}; // center of rotation in airfoil coordinate system

#define angular_rot (t < pitchStartTime? 0: \
                     rotationRate * (1 - exp(-4.6 * (t - pitchStartTime))))

#define pitch_angle (t < pitchStartTime? theta_p: \
                     angular_rot * (t - pitchStartTime) + theta_p)

#define airfoil_thickness(x) (5 * MAX_THICKNESS * ((0.2969*(sqrt(x)))	\
						   -(0.1260*x)		\
						   -(0.3516*(sq(x)))	\
						   +(0.2843*(cube(x)))	\
						   -(0.1036*(pow(x, 4.)))))

face vector muv[];

u.n[left] = dirichlet ((U0));
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);


u_x_ibm_dirichlet (angular_rot * (y - ci.y))
u_y_ibm_dirichlet (-angular_rot * (x - ci.x))


void airfoil_shape (scalar c, face vector f, double theta, vertex scalar phii = {0})
{

  vertex scalar phi = automatic (phii);
  
  double chord = CHORD_LENGTH;
  
  foreach_vertex() {
    double XX = (cr.x + (x - ci.x)*cos (theta) - (y - ci.y)*sin (theta)) / chord;
    double YY = (cr.y + (x - ci.x)*sin (theta) + (y - ci.y)*cos (theta)) / chord;
    
    if (XX < 0) { // leading edge cap
      phi[] = sq(XX - rr) + sq(YY) - sq(rr);
    }
    else if (XX >= 0. && XX < trailingOffset) {
      double yt = airfoil_thickness(XX);
      phi[] = sq(YY) - sq(yt);
    }
    else if (XX >= trailingOffset && XX <= chord) { // trailing edge cap
      double trailingRadius = airfoil_thickness(trailingOffset) + 0.0000171016;
      phi[] = sq(XX - trailingOffset) + sq(YY) - sq(trailingRadius);
    }
    else {
      phi[] = 1.;
    }
  } 
  boundary ({phi});
  fractions (phi, c, f);
}


int main(){
  size(L0);
  init_grid (1 << (MIN_LEVEL));
  mu = muv;
  TOLERANCE = 1.e-5 [*];

  run();
}


event init (t = 0) {
  // initial mesh refinement
  scalar ibmsf[];
  astats ss;
  int ic = 0;
  do {
    ic++;
    airfoil_shape (ibm, ibmf, theta_p);
    foreach()
      ibmsf[] = vertex_average(point, ibm);
    ss = adapt_wavelet ({ibmsf}, (double[]) {1.e-30},
			maxlevel = LEVEL, minlevel = MIN_LEVEL);
  } while ((ss.nf || ss.nc) && ic < 100);
  
  airfoil_shape(ibm, ibmf, theta_p);
  foreach() {
    u.x[] = ibm[]*U0;
  }
  boundary((scalar *){u});
}

event moving_airfoil (i++)
{
  double pitchAngle = pitch_angle;
  airfoil_shape(ibm, ibmf, pitchAngle);   
}


event properties (i++) {
  foreach_face() {
    muv.x[] = fm.x[]*(U0)*(CHORD_LENGTH)/(Re);
  }
  boundary ((scalar *) {muv});
}

scalar pid[];
event logfile (i++; t <= t_end) {

  coord Fp, Fmu;
  ibm_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(CHORD_LENGTH));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(CHORD_LENGTH));

  foreach()
    pid[] = pid();
 
  double pitchAngle = pitch_angle * 180 / M_PI;
  fprintf (stderr, "%d %g %d %d %d %d %g %g %g\n",
	   i, t, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, pitchAngle);
  
}

event screenshot (t = {0, 1.4831, 1.8}) {
  scalar omega[];
  vorticity (u , omega);

  char name[80];
  sprintf (name, "pres-%g.png", t);
  FILE * fp1 = fopen (name, "w");
  view (fov = 2, camera = "front", 
        tx = -0.25, ty = -0.5, bg = {1,1,1},
	width = 3000, height = 3000);
  clear();
  draw_vof ("ibm", "ibmf", filled = -1, lw = 3);
  squares ("omega", map = cool_warm, min = -100, max = 100);
  save (fp = fp1);
  fclose (fp1);
}

event adapt (i++) {
  scalar ibmsf[];
  foreach()
    ibmsf[] = vertex_average(point, ibm);
  adapt_wavelet ({ibmsf,u}, (double[]){1.e-15,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = MIN_LEVEL);
  event("vof");
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%g\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fprintf (fp, "%d\t%d\t%d\t%d\n", mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  fflush (fp);
  return 1;
}
