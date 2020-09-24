//============================================================================
// Name        : external_force.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include<math.h>

#include "param.inc"

void ext_force(	const int ntotal,
				const double mass[],
				const double x[][],
				const int niac,
				const int pair_i[],
				const int pair__j[],
				const int itype[],
				const double hsml,
				double dvxdt[][] ){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the external forces, e.g. gravitational forces,
//	The forces from the interactions with boundary virtual particles
//	are also calculated here as external forces.

// 	ntotal 	: Number of particles 							[in]
// 	mass	: Particle masses 								[in]
// 	x		: Coordinates of all particles					[in]
//	pair_i 	: List of first partner of interaction pair		[in]
// 	pair_j 	: List of second partner of interaction pair	[in]
// 	itype	: type of particles								[in]
// 	hsml	: Smoothing Length								[in]
//	dvxdt	: Acceleration with respect to x, y and z		[out]
//----------------------------------------------------------------------------

	double dx(dim), rr, f, rrO, dd, p1, p2;

	for( int i = 0; i < ntotal; i++ ){
		for( int d = 0; d < dim; d++ ){
			dvxdt(d, i) = .0;
		}
	}

	// Consider self-gravity or not ?
	if (self_gravity){
		for( int i = 0; i < ntotal; i++){
			dvxdt(dim, i) = -9.8;
		}
	}
	// Boundary particle force and penalty anti-penetration force.
	rrO = 1.25e-5;
	dd 	= 1.0e-2;
	p1 	= 12.0;
	p2 	= 4.0;

	for( int k = 0; k < niac; k++ ){

		const int i = pair_i(k);
		const int j = pair_j(k);

		if( itype(i) > 0 && itype(j) < O ){
			rr = .0;

			for( int d = 0; d < dim; d++){
				dx(d) 	= x(d,i) - x(d,j);
				rr 		= rr + dx(d)*dx(d);
			}

			rr = sqrt(rr);

			if(rr < rrO){
				f = ( pow(rrO/rr,p1) - pow(rrO/rr,p2))/(rr*rr);

				for( int d = 0; d < dim; d++){
					dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f;
				}
			}

		}
	}

}











