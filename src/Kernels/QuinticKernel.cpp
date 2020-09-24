/*
 * QuinticKernel.cpp
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */

#include "QuinticKernel.h"

namespace SPH {

template < uint8_t DIM >
QuinticKernel<DIM>::~QuinticKernel() {
	// TODO Auto-generated destructor stub
}

template < uint8_t DIM >
QuinticKernel<DIM>::QuinticKernel(const QuinticKernel &other) {
	// TODO Auto-generated constructor stub

}

template < uint8_t DIM >
QuinticKernel<DIM>::QuinticKernel(QuinticKernel &&other) {
	// TODO Auto-generated constructor stub

}

template < uint8_t DIM >
QuinticKernel<DIM>& QuinticKernel<DIM>::operator=(const QuinticKernel &other) {
	// TODO Auto-generated method stub

}

template < uint8_t DIM >
QuinticKernel<DIM>& QuinticKernel<DIM>::operator=(QuinticKernel &&other) {
	// TODO Auto-generated method stub

}

template < uint8_t DIM >
constexpr double QuinticKernel<DIM>::k() const {
	return 3.0;
}

template < uint8_t DIM >
double QuinticKernel<DIM>::W() const {
	const double f 		= factor();
	constexpr double q	= .0;

	return f * ( pow( 3.0 - q , 5.0 ) - 6.0 * pow( 2.0 - q , 5.0 ) + 15.0 * pow( 1.0 - q , 5.0 ) );

}

template < uint8_t DIM >
double QuinticKernel<DIM>::W( const double r ) const {

	double W = .0;
	const double f 	= factor();
	const double q	= r / this->hsml();

	if( q >= .0 && q <= 1.0 ){
		W = f * ( pow( 3.0 - q , 5.0 ) - 6.0 * pow( 2.0 - q , 5.0 ) + 15.0 * pow( 1.0 - q , 5.0 ) );
	}
	else if( q > 1.0 && q <= 2.0 ){
		W = f * ( pow( 3.0 - q , 5.0 ) - 6.0 * pow( 2.0 - q , 5.0 ) );
	}
	else if( q > 2.0 && q <= 3.0 ){
		W = f * pow( 3.0 - q , 5 );
	}

	return W;
}

template < uint8_t DIM >
constexpr vDdf<DIM> QuinticKernel<DIM>::dWdx() const {
	return vDdf<DIM>();
}

template < uint8_t DIM >
vDdf<DIM> QuinticKernel<DIM>::dWdx( const double r , const vDdf<DIM> & dx ) const {

	const double q	= r / this->hsml();

	double factor = .0;
	if( q > .0 && q <= 1.0 ){
		factor = this->f() * ( ( -120.0 + 120.0*q - 50.0 * pow( q , 2.0 ) ) / pow(  this->hsml() , 2.0 ));
	}
	else if( q > 1.0 && q <= 2.0 ){
		factor = this->f() * ( -5.0 * pow( 3.0 - q , 4.0 ) + 30.0 * pow( 2.0 - q , 4.0 ) )/  this->hsml() *  r;
	}
	else if( q > 2.0 && q <= 3.0 ){
		factor = this->f() * ( -5.0 * pow( 3.0 - q , 4.0 ) ) /  this->hsml() * r;
	}

	return factor * dx;
}

template <>
double QuinticKernel<1>::factor () const {
	return 1.0 / ( 120.0 * hsml() );
}

template <>
double QuinticKernel<2>::factor () const {
	return 7.0 / ( 478.0 * CONST::PI * pow( hsml() , 2.0 ) );
}

template <>
double QuinticKernel<3>::factor () const {
	return 1.0 / ( 120.0 * CONST::PI * pow( hsml() , 3.0 ) );
}


template class QuinticKernel<1>;
template class QuinticKernel<2>;
template class QuinticKernel<3>;

} /* namespace SPH */

































