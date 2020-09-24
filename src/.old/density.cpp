//============================================================================
// Name        : density.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include "param.inc"

void sum_density(	const int ntotal,
					const double hsml[],
					const double mass[],
					const int niac,
					const int pair_i[],
					const int pair_j[],
					const double w[],
					const int itype[],
					double rho[]
							   ){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the density with SPH summation algorithm.

// 	ntotal 	: Number of particles 							[in]
// 	hsml	: Smoothing Length 								[in]
// 	mass 	: Particle masses 								[in]
// 	niac 	: Number of interaction pairs 					[in]
// 	pair_i 	: List of first partner of interaction pair 	[in]
// 	pair_j 	: List of second partner of interaction pair 	[in]
// 	w 		: Kernel for all interaction pairs 				[in]
// 	itype 	: type of particles 							[in]
// 	x 		: Coordinates of all particles 					[in]
// 	rho 	: Density 										[out]
//----------------------------------------------------------------------------

	double selfdens, hv(dim), r, wi(maxn);


	// wi(maxn) --- integration'of the kernel itself
	for( int d = 0; d < dim; d++){
		hv(d) = .0;
	}

	// Self density of each particle: Wii (Kernel for distance 0) and take contribution of particle itself:
	r = .0;

	// Firstly calculate the integration of the kernel over the space
	for( int i = 0; i < ntotal; i++){
		kernel(r, hv, hsml(i), selfdens, hv);
		wi(i) 	= selfdens*mass(i)/rho(i);
	}

	for( int k = 0; k < niac; k++){
		const int i = pair_i(k);
		const int j = pair_j(k);
		wi(i) = wi(i) + mass(j)/rho(j)*w(k);
		wi(j) = wi(j) + mass(i)/rho(i)*w(k);
	}

	// Secondly calculate the rho integration over the space
	for( int i = 0; i < ntotal; i++){
		kernel(r, hv, hsml(i), selfdens, hv);
		rho(i) = selfdens*mass(i);
	}

	// Calculate SPH sum for rho
	for( int k = 0; k < niac; k++){
		const int i = pair_i(k);
		const int j = pair_j(k);
		rho(i) = rho(i) + mass(j)*w(k);
		rho(j) = rho(j) + mass(i)*w(k);
	}

	// Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
	if (nor_density){

		for(int i = 0; i < ntotal; i++){
			rho(i) = rho(i)/wi(i);
		}
	}

}


void con_density(	const int ntotal,
					const double mass[],
					const int niac,
					const int pair_i[],
					const int pair_j[],
					const double dwdx[][],
					const double vx[][],
					const int itype[],
					const double x[][],
					const double rho[],
					double drhodt[]){

//----------------------------------------------------------------------------
// 	Subroutine to calculate 'the density with SPH continuity approach.

// 	ntotal	: Number of particles 								[in]
// 	mass	: Particle masses 									[in]
// 	niac	: Number of interaction pairs 						[in]
// 	pair_i	: List of first partner of interaction pair 		[in]
// 	pair_j	: List of second partner of interaction pair		[in]
// 	dwdx	: derivation of Kernel for all interaction pairs 	[in]
// 	vx		: Velocities of all particles 						[in]
// 	itype	: type of particles 								[in]
// 	x		: Coordinates of all particles 						[in]
// 	rho		: Density 											[in]
// 	drhodt	: Density change rate of each particle 				[out]
//----------------------------------------------------------------------------

	double vcc, dvx(dim);

	for( int i = 0; i < ntotal; i++){
		drhodt(i) = .0;
	}

	for( int k = 0; k < niac; k++){

		const int i = pair_i(k);
		const int j = pair_j(k);

		for( int d = 0; d < dim; d++){
			dvx(d) = vx(d,i) - vx(d,j);
		}

		vcc = .0;
		for( int d = 0; d < dim; d++){
			vcc = vcc + dvx(d)*dwdx(d,k);
		}

		drhodt(i) = drhodt(i) + mass(j)*vcc;
		drhodt(j) = drhodt(j) + mass(i)*vcc;
	}

}














