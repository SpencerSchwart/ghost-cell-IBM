# General Information
This repo implements a ghost cell immersed boundary method (IBM) for the opensource computational fluid dynamics (CFD) software [Basilisk](http://basilisk.fr).

# Method
## Cell Metrics
The most complicated part of this method is the various cell metrics it introduces. In [Basilisk](http://basilisk.fr), cell center metric, *cm*, and face metric, **fm**, are primarily used to scale the solution to account for an axisymmetric grid or other non-uniform meshes. In normal 2D/3D cartesian coordinates, fm = cm = 1. When used with Basilisk's cut cell method, [embedded boundaries](http://basilisk.fr/src/embed.h), the cell metrics are set to equal the liquid volume fraction fields

```c
fm = fs;
cm = cs;
```
where fs = cs = 1 if cell is full of *liquid* and fs = cs = 0 for pure *solid* cells. Likewise, 0 < cs or fs < 1 for fractional, interfacial cells.

Although the cut cell method is relatively complex, it benefits from consistent use of fm and cm. However, our ghost cell IBM uses methodology from both standard ghost cell and cut cell methods. This hybrid approach complicates the otherwise trivial usage of cell metrics.

The present ghost cell IBM uses a standard ghost cell approach to solve the advection and viscous terms in the Navier-Stokes equations and uses the cut cell method to solve the pressure-Poisson equation and in the two-phase solver. This means, in the code, we make changes to several Basilisk's fields, as listed below:

```c
fm = ibmFaces;
cm = ibmCells;

alpha = ibmf / density; // density = 1 for single phase
rho = ibm * density;
```
where ibm and ibmf are analgous to cs and fs, respectively. **ibmFaces** and *ibmCells* are defined as follows

```c
ibmFaces = {
    1, if cell is a ghost cell or is mostly/all liquid (ibm > 0.5)
    0, otherwise
}

ibmCells = {
    1, if cell is mostly/all liquid (ibm > 0.5)
    0, otherwise
}
```
Considering the above definitions, fm and cm are used primarily to solve the advection and diffusion terms while alpha and rho are primarily used to solve the pressure-Poisson equation and the two-phase solver. This way, we avoid solving the intermediate velocity field, **u<sup>*</sup>**, using the cut cell method, which adds many complexities such as special treatment of small cells and limiting the timestep size. Moreover, we ensure we are solving the pressure-Poisson equation using cut cell techniques, which has been shown to reduce spurious, unphysical pressure oscillations.
