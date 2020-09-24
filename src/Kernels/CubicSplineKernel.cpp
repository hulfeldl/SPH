/*
 * CubicSplineKernel.cpp
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */




#include "CubicSplineKernel.h"

namespace SPH {

template < uint8_t DIM >
CubicSplineKernel<DIM>::~CubicSplineKernel() {
	// TODO Auto-generated destructor stub
}

template < uint8_t DIM >
CubicSplineKernel<DIM>::CubicSplineKernel(CubicSplineKernel &&other) {
	// TODO Auto-generated constructor stub

}

template < uint8_t DIM >
CubicSplineKernel<DIM>& CubicSplineKernel<DIM>::operator=( const CubicSplineKernel &other) {
	return *this = other;

}

template < uint8_t DIM >
CubicSplineKernel<DIM>& CubicSplineKernel<DIM>::operator=(CubicSplineKernel &&other) {
	// TODO Auto-generated method stub

}

template < uint8_t DIM >
void CubicSplineKernel<DIM>::update_f() {
	this->set_f( factor() );
}

template < uint8_t DIM >
double CubicSplineKernel<DIM>::W() const {

	return this->f() * ( 2.0 / 3.0 );
}

template < uint8_t DIM >
double CubicSplineKernel<DIM>::W( const double r  ) const {

	const double q 	= r / this->hsml();

	double W = .0;
	if ( q >= 0.0 && q <= 1.0 ){
		W = this->f() * ( 2.0 / 3.0 - pow( q , 2.0 ) + 0.5 * pow( q , 3.0 ) );
	}
	else if ( q > 1.0 && q <= 2.0 ){
		W = this->f() * 1.0 / 6.0 * pow( 2.0 - q , 3.0 );

	}

	return W;
}

template < uint8_t DIM >
vDdf<DIM> CubicSplineKernel<DIM>::dWdx( const double r , const vDdf<DIM> & dx ) const {

	const double q 	= r / this->hsml();

	vDdf<DIM> dWdx;
	double factor = .0;
	if ( q > 0.0 && q <= 1.0 ){
		factor = this->f() * ( -2.0 + 3.0/2.0*q ) / pow( this->hsml() , 2.0 );
	}
	else if ( q > 1.0 && q <= 2.0 ){
		factor = - this->f() * 1.0 / 6.0 * 3.0 * pow( 2.0 - q , 2.0 ) / ( this->hsml() * r );
	}

	return factor * dx;
}

template < uint8_t DIM >
double CubicSplineKernel<DIM>::factor () const {
	MYASSERT (false, "Dimension must be 1, 2 or 3!");
	return .0;
}

template <>
double CubicSplineKernel<1>::factor () const {
	return 1.0 / this->hsml();
}

template <>
double CubicSplineKernel<2>::factor () const {
	return 15.0 /( 7.0 * CONST::PI * pow( this->hsml() , 2.0 ) );
}

template <>
double CubicSplineKernel<3>::factor () const {
	return 3.0 / ( 2.0 * CONST::PI * pow( this->hsml() , 3.0 ) );
}

template class CubicSplineKernel<1>;
template class CubicSplineKernel<2>;
template class CubicSplineKernel<3>;

} /* namespace SPH */










