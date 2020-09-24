//============================================================================
// Name        : art_heat.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================

#include "param.inc"

void art_heat(	const int ntotal,
				const double hsml[],
				const double mass[],
				const double x[][],
				const double vx[][],
				const int niac,
				const double rho[],
				const double u[],
				const double c[],
				const int pair_i[],
				const int pair_j[],
				const double w[],
				const double dwdx[],
				double dedt[]){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the artificial heatl Fulk, 1994, p, a-17)

// 	ntotal	: Number of particles								[in]
// 	hsml	: Smoothing Length									[in]
// 	mass	: Particle masses									[in]
// 	x		: Coordinates of all particles						[in]
// 	vx		: Velocities of all particles						[in]
// 	rho		: Density											[in]
// 	u		: specific internal energy							[in]
// 	c		: Sound velocity									[in]
// 	niac	: Number of interaction pairs						[in]
//	pair_i	: List of first partner of interaction pair			[in]
//	pair_j	: List of second partner of interaction pair		[in]
//	w		: Kernel for all interaction pairs					[in]
//	dwdx	: Derivative of kernel with respect to x, y and z	[in]
//	dedt	: produced artificial heat, adding to energy Eq.	[out]

//----------------------------------------------------------------------------


	double dx, dvx(dim), vr, rr, h, mc, mrho, mhsml, vcc(maxn), hvcc, mui, muj, muij, rdwdx, g1,g2;


	//------------------------------------------------------------------------
	// Parameter for the artificial heat conduction:
	g1	= 0.1;
	g2	= 1.0;

	for(int i = 0; i < ntotal; i++){
		vcc(i) 	= .0;
		dedt(i) = .0;
	}

	for (int k = 0; k < niac; k++){
		const int i = pair_i(k);
		const int j = pair_j(k);

		for(int d = 0; d < dim; d++){
			dvx(d) = vx(d,j) - vx(d,i);
		}

		hvcc = .0;
		for(int d = 0; d < dim; d++){
			hvcc = hvcc + dvx(d)*dwdx(d,k);
		}

		vcc(i) = vcc(i) + mass(j)*hvcc/rho(j);
		vcc(j) = vcc(j) + mass(i)*hvcc/rho(i);
	}

	for(int k = 0; k < niac; k++){

		const int i = pair_i (k);
		const int j = pair_j(k);

		mhsml	= 0.5*(hsml(i)+hsml(j));
		mrho 	= 0.5*(rho(i) + rho(j));

		rr 		= .0;
		rdwdx 	= .0;
		for(int d = 0; d < dim; d++){
			dx = x(d,i) - x(d,j);
			rr = rr + dx*dx;
			rdwdx = rdwdx + dx*dwdx(d,k);
		}

		mui		= g1*hsml(i)*o(i) + g2*(hsml(i)*hsml(i))*(abs(vcc(i))-vcc(i));
		muj		= g1*hsml(j)*c(j) + g2*(hsml(j)*hsml(j))*(abs(vcc(j))-vcc(j));
		muij	= 0.5*(mui+muj);

		h 		= muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx;
		dedt(i) = dedt(i) + mass(j)*h*(u(i)-u(j));
		dedt(j) = dedt(j) + mass (i) *h* (u (j ) -u(i) );
	}

	for(int i = 0; i < ntotal; i++){
		dedt(i) = 2.0*dedt(i);
	}

}















