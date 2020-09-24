//============================================================================
// Name        : hsml.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <math.h>

#include "param.inc"

void h_upgrade( const double dt,
				const int ntotal,
				const double mass[],
				const double vx[][],
				const double rho[],
				const int niac,
				const int pair_i[],
				const int pair_j[],
				const double dwdx[][],
				double hsml[] ){

//----------------------------------------------------------------------------
//	Subroutine to evolve smoothing length

// 	dt		: time step											[in]
// 	ntotal 	: Number of particles 								[in]
//	mass	: Particle masses									[in]
//	vx		: Velocities of all particles						[in]
//	rho		: Density											[in]
//	niac	: Number of interaction pairs						[in]
//	pair_i	: List of first partner of interaction pair			[in]
//	pair_j 	: List of second partner of interaction pair		[in]
//	dwdx	: Derivative of kernel with respect to x, y and z	[in]
//	hsml	: Smoothing Length									[in/out]
//----------------------------------------------------------------------------

	double fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn);

	if (sle == 0 ){
		// Keep smoothing length unchanged.
		return;
	}
	else if(sle == 1){
		fac = 2.0;
		for(int i = 0; i < ntotal; i++){
			hsml(i) = fac*(pow( (mass(i)/rho(i)), (1.0/double(dim)) ));
		}
	}
	else if (sle == 2){
		// dh/dt = (-I/dim)*(h/rho)*(drho/dt).

		for( int i = 0; i < ntotal; i++){
			vcc(i) = .0;
		}

		for( int k = 0; k < niac; k++){
			const int i = pair_i(k);
			const int j = pair_j(k);

			for( int d = 0; d < dim; d++){
				dvx(d) = vx(d,j) - vx(d,i);
			}

			hvcc = .0;

			for( int d = 0; d < dim; d++){
				hvcc = hvcc + dvx(d)*dwdx(d,k);
			}

			vcc(i) = vcc(i) + mass(j)*hvcc/rho(j);
			vcc(j) = vcc(j) + mass(i)*hvcc/rho(i);
		}

		for( int i = 0; i < ntotal; i++){
			dhsml(i) 	= (hsml(i)/dim)*vcc(i);
			hsml(i)	 	= hsml(i) + dt*dhsml(i);

			if (hsml(i) <= 0){
				hsml(i) = hsml(i) - dt*dhsml(i);
			}
		}
	}

}



