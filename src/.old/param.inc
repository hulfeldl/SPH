//============================================================================
// Name        : param.inc
// Author      : Lorenz
// Version     :
// Copyright   : Your copyright notice
// Description : param.inc from "10.4.5 List of the FORTRAN source files"
//============================================================================


constexpr double pi = 3.14159265358979323846;

//----------------------------------------------------------------------------
//	Including file for parameters and constants used
// 	in the entire SPH software packages.
//----------------------------------------------------------------------------

// dim: Dimension of the problem (1, 2 or 3)
int dim = 1;

// maxn: Maximum number of particles
// max_interation : Maximum number of interaction pairs
int maxn = 12e2;
int max_interaction = 100*maxn;

// Parameters for the computational geometry,
// x_maxgeom: Upper limit of allowed x-regime
// x_mingeom: Lower limit of allowed x-regime
// y_maxgeom: Upper limit of allowed y-regime
// y_mingeom: Lower limit of allowed y-regime
// z_maxgeom: Upper limit of allowed z-regime
// z_mingeom: Lower limit of allowed z-regime

double x_maxgeom = 10.0, x_mingeom = -10.0;
double y_maxgeom = 10.0, y_mingeom = -10.0;
double z_maxgeom = 10.0, z_mingeom = -10.0;

// SPH algorithm for particle approximation (pa_sph)
// pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
// 			2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
int pa_sph = 2;

// Nearest neighboring particle searching (nnps) method
// nnps =	1 : Simplest and direct searching
// 			2 : Sorting grid linked list
// 			3 : Tree algorithm
int nnps = 1;

// Smoothing length evolution (sle) algorithm
// sle =	0: Keep unchanged
//			1: h 		= fac * (m/rho)*(1/dim)
//			2: dh/dt 	= (-1/dim)*(h/rho)*(drho/dt)
//			3: Other approaches (e.g. h = h_0 • (rho_0/rho)**(I/dim) )
int sle = 0;

// Smoothing kernel function
// skf = 	1: Cubic spline kernel by W4 - Spline (Monaghan 1985)
//			2: Gauss kernel (Gingold and Monaghan 1981)
// 			3: Quintic kernel (Morris 1997)

int skf = 1;


// Switches for different senarios
//----------------------------------------------------------------------------
// summation_density = 	true: 	Use density summation model in the code
// 						false: 	Use continuiity equation
// average_velocity = 	true: 	Monaghan treatment on average velocity,
//						false: 	No average treatment.
// config_input = 		true: 	Load initial configuration data,
// 						false: 	Generate initial configuration.
// virtual_part = 		true: 	Use vritual particle,
// 						false: 	No use of vritual particle.
// vp_input = 			true: 	Load virtual particle information,
// 						false: 	Generate virtual particle information.
// visc = 				true: 	Consider viscosity,
// 						false: 	No viscosity.
// ex_force = 			true: 	Consider external force,
// 						false: 	No external force.
// visc_artificial = 	true: 	Consider artificial viscosity,
// 						false: 	No considering of artificial viscosity.
// heat_artificial = 	true: 	Consider artificial heating,
// 						false: 	No considering of artificial heating.
// self_gravity = 		true: 	Considering self_gravity,
// 						false: 	No considering of self_gravity
// nor_density = 		true: 	Density normalization by using CSPM,
// 						false: 	No normalization.

bool summation_density 	= true;
bool average_velocity 	= false;
bool config_input 		= false;
bool virtual_part 		= false;
bool vp_input 			= false;
bool visc 				= false;
bool ex_force 			= false;
bool visc_artificial 	= true;
bool heat_artificial 	= false;
bool self_gravity 		= false;
bool nor_density 		= false;


// Symmetry of the problem
// nsym =	0: no symmetry
// 			1: axis symmetry
// 			2: center symmetry

int nsym = 0;

// Control parameters for output
// int_sta = true 	Print statistics about SPH particle interactions including virtual particle information.
// print_step: 		Print Timestep (On Screen)
// save_step : 		Save Timestep (To Disk File)
// moni_particle: 	The particle number for information monitoring.

bool int_stat 		= true;
int print_step 		= 100;
int save_step 		= 500;
int moni_particle 	= 1600;

// Simulation cases
// shocktube 	= true : carry out shock tube simulation
// shearcavity 	= true : carry out shear cavity simulation

bool shocktube 		= true;
bool shearcavity 	= false;



















