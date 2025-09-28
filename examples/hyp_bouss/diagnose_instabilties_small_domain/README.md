## Minimal example to understand the source of instabilities when solving the hyperbolic approximation of SGN
## by Escalante et al. (2019).


##### Description:

This example only considers a small domain around the Galapagos Islands. 
Some instabilties are observed when the relaxation parameter ($c^2$) increases.
As $c^2\to\infty$ we should recover SGN, and as $c^2\to 0$ we should recover the shallow water equations.

Some things we know:
- The issue is not related to AMR. It seems like it was triggered more easily before with the inconsistency at ghost cells between patches, but  can be observed even with a single grid if $c^2$ is taken large enough.
- The instabilities seem to be introduced by the convective terms (the normal Riemann solver). Additional diffusion, like a smaller CFL number or using an approximate integrator for the source terms (e.g. backward Euler instead of the exact one) seems to slightly help, but not too much.
- The instabilities arise also in Cartesian coordinates and with no transition to SWE (i.e. constant but large $c^2$), the determining factor seems to be large jumps in the bathymetry (for instance, one can use the artificial bathymetry etopo_no_dry.asc such that there is never a transition). In particular, this points to a non-conservative term appearing in the last equation of the system, proportional to $c^2 b_x u$. If we artificially remove this term from the equations, the instabilities seem to dissapear (see the local bouss_module.f90 and the one in src/2d/hyp_bouss). Moreover, in different formulations of the hyperbolic approximation (e.g. Duran & Richard 2024) this term is not present in the equations. Could it be possible that, in non-dimensional variables and before taking a hyperbolic approximation, this term is proportional to (shallow regime) terms that were neglected in the derivation of the model, and the hyperbolic appoximation amplifies it unphysically?


##### To run this example:
make data
make .output
make plots


##### Local source code:
The local source files are only modified to rule out certain possibilities. For instance, in this example, the aux array  which indicates what system we are solving and the local value (can be made space dependent) of $c^2$ is fixed in time (the code in src/2d/hyp_bouss sets it dynamically and avoids recomputation like its done for the bahymetry). Also, the non-conservative term mentioned above is commented out in the local bouss_module.f90 (which contains the Riemann solvers).