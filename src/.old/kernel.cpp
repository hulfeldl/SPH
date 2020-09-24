//============================================================================
// Name        : kernel.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>
#include <math.h>

#include "param.inc"

void kernel(const double r,
			const double dx[],
			const double hsml,
			double w,
			double dwdx[] ){

//----------------------------------------------------------------------------
// 	Subroutine to calculate the smoothing kernel wij and its
// 	derivatives dwdxij.
// 	if skf 	= 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
// 			= 2, Gauss kernel (Gingold and Monaghan 1981)
// 			= 3, Quintic kernel (Morris 1997)

// 	r 		: Distance between particles i and j				[in]
// 	dx		: x-, y- and z-distance between i and j 			[in]
// 	hsml	: Smoothing length									[in]
//	w		: Kernel for all interaction pairs					[out]
// 	dwdx	: Derivative of kernel with respect to x, y and z	[out]
//----------------------------------------------------------------------------

	double q, dw, factor;

	double q = r/hsml;
	w = .0;

	for( int d = 0; d < dim; d++ ){
		dwdx(d) = .0;
	}

	if ( skf == 1){

		if ( dim == 1){
			factor = 1.0/hsml;
		}
		else if ( dim == 2 ){
			factor = 15.0/(7.0*pi*hsml*hsml);
		}
		else if ( dim == 3 ){
			factor = 3.0/(2.0*pi*hsml*hsml*hsml);
		}
		else{
			std::cout << " >>> Error <<< : Wrong dimension: Dim = " << dim;
			return;
		}

		if ( q >= O && q <= 1.0 ){

			w = factor * ( 2.0/3.0 - q*q + 0.5*pow(q,3) );

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = factor * (-2.0 + 3.0/2.0*q)/pow(hsml,2) * dx(d);
			}
		}
		else if ( q > 1.0 && q <= 2.0 ){

			w = factor * 1.0/6.0 * pow( 2.0 - q , 3 );

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = -factor * 1.0/6.0 * 3.0 * pow( 2.0 - q , 2)/hsml * (dx(d)/r);
			}

		}
		else{

			w = .0;

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = .0;
			}
		}
	}
	else if ( skf == 2 ){

		factor = 1.0 / ( pow(hsml,dim) * pow( pi , dim/2.0 ) );

		if ( q >= .0 && q <= 3.0 ){

			w = factor * exp(-q*q);

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = w * ( - 2.0 * dx(d)/hsml/hsml);
			}
		}
		else{
			w = .0;

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = .0;
			}
		}

	}
	else if ( skf == 3 ){

		if ( dim == 1){
			factor = 1.0 / ( 120.0*hsml );
		}
		else if ( dim == 2 ){
			factor = 7.e0 / ( 478.0*pi*hsml*hsml );
		}
		else if ( dim == 3 ){
			factor = 1.0 / ( 120.0*pi*hsml*hsml*hsml );
		}
		else{
			std::cout << " >>> Error <<< : Wrong dimension: Dim = " << dim;
			return;
		}

		if( q >= .0 && q <= 1.0 ){

			w = factor * ( pow( 3.0 - q , 5 ) - 6.0*pow( 2.0 - q , 5 ) + 15.0*pow( 1.0 - q , 5 ) );

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = factor * ( ( -120.0 + 120.0*q - 50.0 * pow(q,2) )/ pow(hsml,2) * dx(d) );
			}
		}
		else if( q > 1.0 && q <= 2.0 ){

			w = factor * ( pow( 3.0 - q , 5 ) - 6.0 * pow( 2.0 - q , 5 ) );

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = factor * ( -5.0 * pow( 3.0 - q , 4 ) + 30.0 * pow( 2.0 - q , 4 ) )/ hsml * ( dx(d) / r );
			}
		}
		else if( q > 2.0 && q <= 3.0 ){

			w = factor * pow( 3.0 - q , 5 );

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = factor * ( -5.0 * pow( 3.0 - q , 4 ) ) / hsml * ( dx(d) / r );
			}
		}
		else{

			w = .0;

			for( int d = 0; d < dim; d++ ){
				dwdx(d) = .0;
			}

		}

	}

}








