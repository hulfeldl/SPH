//============================================================================
// Name        : av_vel.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include "param.inc"

void av_vel(const int ntotal,
			const double mass[],
			const int niac,
			const int pair_i[],
			const int pair_j[],
			const double w[],
			const double vx[][],
			const double rho[],
			double av[][]){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the average velocity to correct velocity
// 	for preventing.penetration (monaghan, 1992)

// 	ntotal	: Number of particles 							[in]
// 	mass	: Particle masses								[in]
// 	niac	: Number of interaction pairs					[in]
// 	pair_i	: List of first partner of interaction pair		[in]
// 	pair_j	: List of second partner of interaction pair	[in]
// 	w		: Kernel for all interaction pairs				[in]
// 	vx		: Velocity of each particle						[in]
// 	rho		: Density of each particle						[in]
// 	av		: Average velocity of each particle				[out]
//----------------------------------------------------------------------------

	double vcc, dvx(dim), epsilon;

	// epsilon --- a small constants chosen by experence, may lead to instability.
	// for example, for the 1 dimensional shock tube problem, the E <= 0.3

	epsilon = 0.3;

	for( int i = 0; i < ntotal; i++){
		for( int d = 0; d < dim; d++){
			av(d,i) = .0;
		}
	}

	for( int k = 0; k < niac; k++){
		const int i = pair_i(k);
		const int j = pair_j(k);

		for( int d = 0; d < dim; d++){
			dvx(d) 	= vx(d,i) - vx(d,j);
			av(d,i) = av(d,i) - 2*mass(j)*dvx(d)/(rho(i)+rho(j))*w(k);
			av(d,j) = av(d,j) + 2*mass(i)*dvx(d)/(rho(i)+rho(j))*w(k);
		}
	}

	for( int i = 0; i < ntotal; i++){
		for( int d = 0; d < dim; d++){
			av(d,i) = epsilon * av(d,i);
		}
	}

}
















