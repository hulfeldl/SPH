//============================================================================
// Name        : art_visc.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================

#include "param.inc"

void art_visc( 	const int ntotal,
				const double hsml[],
				const double mass[],
				const double x[][],
				const double vx[][],
				const int niac,
				const double rho[],
				const double c[],
				const int pair_i[],
				const int pair_j[],
				const double w[],
				const double dwdx[][],
				double dvxdt[][],
				double dedt[]){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the artificial viscosity (Monaghan, 1992)

// 	ntotal	: Number of particles (including virtual particles)	[in]
// 	hsml	: Smoothing Length									[in]
// 	mass	: Particle masses									[in]
//	x		: Coordinates of all particles						[in]
//	vx		: Velocities of all particles						[in]
// 	niac	: Number of interaction pairs						[in]
// 	rho		: Density											[in]
// 	c		: Temperature										[in]
// 	pair_i	: List of first partner of interaction pair			[in]
// 	pair_j	: List of second partner of interaction pair		[in]
// 	w		: Kernel for all interaction pairs					[in]
// 	dwdx	: Derivative of kernel with respect to x, y and z	[in]
// 	dvxdt	: Acceleration with respect to x, y and z			[out]
// 	dedt	: Change of specific internal energy				[out]
//----------------------------------------------------------------------------

	double dx, dvx(dim), alpha, beta, etq, piv,muv, vr, rr, h, mc, mrho, mhsml;

	// Parameter for the artificial viscosity:
	// Shear viscosity

	alpha 	= 1.0;

	// Bulk viscosity
	beta 	= 1.0;

	// Parameter to avoid singularities
	etq		= .1;

	for(int i = 0; i < ntotal; i++){
		for(int d = 0; d < dim; d++){
			dvxdt(d,i) = .0;
		}

		dedt(i) = .0;
	}


	// Calculate SPH sum for artificial viscosity
	for(int k = 0; k < niac; k++){
		int i = pair_i(k);
		int j = pair_j(k);
		mhsml = 0.5*(hsml(i)+hsml(j));

		vr = .0;
		rr = .0;
		for(int d = 0; d < dim; d++){
			dvx(d) 	= vx(d,i) - vx(d,j);
			dx 		= x(d,i) - x(d,j);
			vr 		= vr + dvx(d)*dx;
			rr 		= rr + dx*dx;
		}


		// Artificial viscous force only if v_ij * r_ij < 0
		if (vr < .0){

			// Calculate muv_ij = hsml v_ij * r_ij / ( r_ij A 2 + hsml*2 etq*2 )
			muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq);

			// Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij*2) / rho_ij
			mc 		= 0.5*(c(i) + c(j));
			mrho 	= 0.5*(rho(i) + rho(j));
			piv 	= (beta*muv - alpha*mc)*muv/mrho;


			// Calculate SPH sum for artificial viscous force
			for(int d = 0; d < dim; d++){
				h 			= -piv*dwdx(d,k);
				dvxdt(d,i) 	= dvxdt(d,i) + mass(j)*h;
				dvxdt(d,j) 	= dvxdt(d,j) - mass(i)*h;
				dedt(i) 	= dedt(i) - mass(j)*dvx(d)*h;
				dedt(j) 	= dedt(j) - mass(i)*dvx(d)*h;
			}
		}
	}


	// Change of specific internal energy
	for(int i = 0; i < ntotal; i++){
		dedt(i) = 0.5*dedt(i);
	}

}










