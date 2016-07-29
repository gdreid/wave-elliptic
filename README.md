# wave-elliptic
Uses Dr. Pretorius' PAMR and AMRD code to solve the wave equation in cylindrical coordinates. This is used as a source for Poisson's equation in Cartesian coordinates. When adaptive mesh refinement is used  or the problem is split up over several processors, the MG variables become non-smooth at the refinement boundaries. It should be noted that the wave equation is not coupled to the potential field; it only sources it.


****

The example in this folder solves the wave equation in cylindrical coordinates:

diff(f,t)   = f_t,
diff(f_t,t) = 1/r diff(r*diff(f,r),r) + diff(f,z$2),

with Neumann BCs on the r=0 boundary and Dirichlet BCs elsewhere. Additionally it solves the elliptic problem:

diff(psi_t0,r$2) + diff(psi_t0,r$2) - f*f = 0

with Neumann BCs on the r=0 boundary and Dirichlet BCs elsewhere.

****

