/*
 * GaussKernel.cpp
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */

#include "GaussKernel.h"

namespace SPH {

template < uint8_t DIM >
GaussKernel<DIM>::~GaussKernel() {
	// TODO Auto-generated destructor stub
}

template < uint8_t DIM >
GaussKernel<DIM>::GaussKernel(const GaussKernel &other) {
	// TODO Auto-generated constructor stub

}

template < uint8_t DIM >
GaussKernel<DIM>::GaussKernel(GaussKernel &&other) {
	// TODO Auto-generated constructor stub

}

template < uint8_t DIM >
GaussKernel<DIM>& GaussKernel<DIM>::operator=(const GaussKernel &other) {
	// TODO Auto-generated method stub

}

template < uint8_t DIM >
GaussKernel<DIM>& GaussKernel<DIM>::operator=(GaussKernel &&other) {
	// TODO Auto-generated method stub

}

template < uint8_t DIM >
constexpr double GaussKernel<DIM>::k() const {
	return 3.0;
}

template < uint8_t DIM >
double GaussKernel<DIM>::W() const {
	return this->f();
}

template < uint8_t DIM >
double GaussKernel<DIM>::W( const double r ) const {
	const double q = r / this->hsml();

	double W = .0;
	if ( q >= .0 && q <= 3.0 ){
		W = this->f() * exp( -q * q );
	}

	return W;
}

template < uint8_t DIM >
constexpr vDdf<DIM> GaussKernel<DIM>::dWdx() const {
	return vDdf<DIM>();
}

template < uint8_t DIM >
vDdf<DIM> GaussKernel<DIM>::dWdx( const double r , const vDdf<DIM> & dx ) const {

	const double q	= r / this->hsml();

	vDdf<DIM> dWdx;
	double factor = .0;
	if ( q > .0 && q <= 3.0 ){
		factor = this->f() * exp( -q * q ) * ( - 2.0  / pow( this->hsml() , 2.0 ) );
	}

	return factor * dx;
}

template < uint8_t DIM >
double GaussKernel<DIM>::factor () const {
	return	1.0 / ( pow( this->hsml(), static_cast<double>( DIM ) ) * pow( CONST::PI , static_cast<double>( DIM ) / 2.0 ) );
}

template class GaussKernel<1>;
template class GaussKernel<2>;
template class GaussKernel<3>;

} /* namespace SPH */



















