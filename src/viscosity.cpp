//============================================================================
// Name        : viscosity.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <stdlib.h>

#include "param.inc"

void viscosity(	const int ntotal,
				const int itype[],
				const double x[][],
				const double rho[],
				double eta[] ){

//----------------------------------------------------------------------------
// 	Subroutine to define the fluid particle viscosity

// 	ntotal	: Number of particles			[in]
//	itype	: Type of particle				[in]
//	x		: Coordinates of all particles	[in]
//	rho		: Density						[in]
//	eta		: Dynamic viscosity				[out]
//----------------------------------------------------------------------------


	for( int i = 0; i < ntotal; i++){

		if ( abs(itype(i)) == 1 ){
			eta(i) = .0;
		}
		else if ( abs(itype(i)) == 2 ){
			eta(i) = 1.0e-3;
		}
	}

}




