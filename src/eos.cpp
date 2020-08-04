//============================================================================
// Name        : eos.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <cmath>
#include <math>

void p_gas(	const double rho,
			const double u,
			double & p,
			double & c ){

//----------------------------------------------------------------------------
// 	Gamma law EOS: subroutine to calculate the pressure and sound

//	rho	: Density 			[in]
//	u	: Internal energy	[in]
//	p	: Pressure			[out]
//	c	: sound velocity	[out]
//----------------------------------------------------------------------------

	// For air (ideal gas)
	const double gamma	= 1.4
	p 		= (gamma - 1)*rho*u;
	c 		= sqrt( (gamma - 1)*u );

}

void p_art_water(	const double rho,
					const double u,
					double & p,
					double & c ){

//----------------------------------------------------------------------------
// 	Artificial equation of state for the artificial compressibility

// 	rho : Density 			[in]
// 	u 	: Internal energy 	[in]
// 	p 	: Pressure 			[out]
// 	c 	: sound velocity 	[out]

// 	Equation of state for artificial compressibility
//----------------------------------------------------------------------------

	double gamma, rho0, b;

	// Artificial EOS, Form 1 (Monaghan, 1994)
	gamma	= 7.0;
	rho0 	= 1000.0;
	b 		= 1.013e5;
	p 		= b*( std::pow(rho/rho0,gamma) - 1);
	c 		= 1480;

	// Artificial EOS, Form 2 (Morris, 1997)
	c  		= 0.01;
	p 		= c*c*rho

}









