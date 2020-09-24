//============================================================================
// Name        : sph.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>

#include "../.old/param.inc"
#include "../.old/time_elapsed.cpp"
#include "../.old/time_integration.cpp"
#include "../.old/time_print.cpp"
#include "input.cpp"
#include "output.cpp"


void SPH(){

//----------------------------------------------------------------------------
// 	This is a three dimensional SPH code, the followings are the
// 	basic parameters needed in this code or calculated by this code

// 	mass	: mass of particles							[in]
//	ntotal	: total particle number						[in]
//	dt		: Time step used in the time integration	[in]
// 	itype	: types of particles						[in]
// 	x		: coordinates of particles					[in/out]
// 	vx		: velocities of particles					[in/out]
//	rho		: dnesities of particles					[in/out]
//	p		: pressure of particles						[in/out]
//	u		: internal energy of particles				[in/out]
//	hsml 	: smoothing lengths of particles			[in/out]
//	c		: sound velocity of particles				[out]
//	s		: entropy of particles						[out]
//	e		: total energy of particles					[out]
//----------------------------------------------------------------------------

	double mass[];
	int ntotal;
	double dt;
	int itype[];
	double x[][];
	double vx[][];
	double rho[];
	double p[];
	double u[];
	double hsml[];
	double c[];
	double s[];
	double e[];

	int maxtimestep, d, m, i, yesorno;

	std::time_t s1, s2;

	time_print();
	time_elapsed(s1);

	if (shocktube){
		dt = 0.005;
	}
	if (shearcavity){
		dt = 5.e-5;
	}

	input(x, vx, mass, rho, p, u, itype, hsml, ntotal);

	yesorno = 1;

	while(yesorno){

		std::cout << "***************************************************" << std::endl;
		std::cout << "        Please input the maximal time steps 		 " << std::endl;
		std::cout << "***************************************************" << std::endl;

		std::cin >> maxtimestep;

		time_integration(x, vx, mass, rho, p, u, c, s, e, itype,hsml, ntotal, maxtimestep, dt );

		output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal);

		std::cout << "***************************************************" << std::endl;
		std::cout << "Are you going to run more time steps? (0=No, 1=yes)" << std::endl;
		std::cout << "***************************************************" << std::endl;

		std::cin >> yesorno;
	}

	time_print();
	time_elapsed(s2);

	std::cout << "Elapsed CPU time = " << s2 - s1 << std::endl;

}

























