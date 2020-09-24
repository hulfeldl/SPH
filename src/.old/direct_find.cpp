//============================================================================
// Name        : direct_find.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <math.h>
#include <iostream>

#include "param.inc"

void direct_find( 	const int itimestep,
					const int ntotal,
					const double hsml[],
					const double x[][],
					int niac,
					int pair_i[],
					int pair_j[],
					double w[],
					double dwdx[][],
					int countiac[] ){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the smoothing funciton for each particle and
// 	the interaction parameters used by the SPH algorithm. Interaction
// 	pairs are determined by directly comparing the particle distance
// 	with the corresponding smoothing length.

// 	itimestep	: Current time step 								[in]
//	ntotal		: Number of particles 								[in]
// 	hsml 		: Smoothing Length 									[in]
//	x			: Coordinates of all particles						[in]
//	niac		: Number of interaction pairs						[out]
//	pair_i		: List of first partner of interaction pair			[out]
//	pair_j		: List of second partner of interaction pair		[out]
//	w			: Kernel for all interaction pairs					[out]
//	dwdx 		: Derivative of kernel with respect to x, y and z 	[out]
//	countiac	: Number of neighboring particles					[out]
//----------------------------------------------------------------------------

	int sumiac, maxiac, miniac, noiac,maxp, minp, scale_k;

	double dxiac(dim), driac, r, mhsml, tdwdx(dim);

	if (skf == 1){
		scale_k = 2;
	}
	else if (skf == 2){
		scale_k = 3;
	}
	else if (skf == 3){
		scale_k = 3;
	}

	for( int i = 0; i < ntotal; i++){
		countiac(i) = 0;
	}

	niac = 0;

	for( int i = 0; i < ntotal-1; i++){

		for( int j = i + 1; j < ntotal; j++){

			driac 		= .0;
			for( int d = 0; d < dim; d++){
				dxiac(d) 	= x(d,i) - x(d,j);
				driac 		= driac + dxiac(d)*dxiac(d);
			}

			mhsml = 0.5*(hsml(i)+hsml(j));

			if (sqrt(driac) < scale_k*mhsml){
				if (niac < max_interaction){

				// Neighboring pair list, and totalinteraction number and
				// the interaction number for each particle
				niac 			= niac + 1;
				pair_i(niac) 	= i;
				pair_j(niac) 	= j;
				r 				= sqrt(driac);
				countiac(i) 	= countiac(i) + 1;
				countiac(j) 	= countiac(j) + 1;


				// Kernel and derivations of kernel
				kernel(r,dxiac,mhsml,w(niac),tdwdx)

				for( int d = 0; d < dim; d++){
					dwdx(d,niac) = tdwdx(d);
				}

				}
				else{
					std::cout << ">>> ERROR <<< : Too many interactions";
					return;
				}

			}
		}
	}

	// Statistics for the interaction
	sumiac 	= 0;
	maxiac 	= 0;
	miniac 	= 1000;
	noiac 	= 0;

	for( int i  = 0; i < ntotal; i++){
		sumiac = sumiac + countiac(i);

		if (countiac(i) > maxiac){
			maxiac = countiac(i);
			maxp = i;
		}

		if (countiac(i) < miniac){
			miniac = countiac(i);
			minp = i;
		}

		if (countiac(i) == 0){
			noiac = noiac + 1;
		}
	}


	if (itimestep % print_step == 0){
		if (int_stat){
			std::cout >> " >> Statistics: interactions per particle:" << endl;
			std::cout >> "**** Particle		: " << maxp << " maximal interactions : " << maxiac << endl;
			std::cout >> "**** Particle		: " << minp << " minimal interactions : " << miniac << endl;
			std::cout >> "**** Average 		: " << real(sumiac)/real(ntotal) << endl;
			std::cout >> "**** Total pairs 	: " << niac << endl;
			std::cout >> "**** Particles with no interactions : " << noiac << endl;
		}
	}

}










