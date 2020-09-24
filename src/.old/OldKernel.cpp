/*
 * Kernel.cpp
 *
 *  Created on: Aug 5, 2020
 *      Author: hulfeldl
 */




#include "OldKernel.h"

#include <math.h>


namespace SPH {

// Abstract Kernel
// ---------------------------------------------------------------------------------------------

template < uint8_t DIM , SmoothingKernel kType >
void AbstractKernel<DIM,kType>::evolveSML( const double drhodt , const double dt ) {

	switch ( m_sle ) {

	case SmoothingLengthEvolution::DensityAlgebraic:
		constexpr double fac 	= 2.0;
		m_Particle.h_densityAlgebraic( fac );
		break;

	case SmoothingLengthEvolution::SmoothingLengthODE:
		m_Particle.h_smlODE( drhodt );
		break;

	case SmoothingLengthEvolution::OtherApproach:

		break;
	}
}

template < uint8_t DIM , SmoothingKernel kType >
void AbstractKernel<DIM,kType>::update_hsml( const double hsml ) {
	m_hsml = hsml;
}

template < uint8_t DIM , SmoothingKernel kType >
void AbstractKernel<DIM,kType>::update_hsml(const double dhsmldt ,
											const double dt ) {
	const double hsml_old = m_hsml;
	m_hsml += dhsmldt * dt;

	if ( m_hsml <= 0 ){
		m_hsml = hsml_old;
	}
}


// Kernel
// ---------------------------------------------------------------------------------------------

template < uint8_t DIM , SmoothingKernel kType >
OldKernel<DIM,kType>::OldKernel() :
	m_sle ( SmoothingLengthEvolution::NONE ),
	m_hsml ( .0 ) {
}

// Std. Constructor
template < uint8_t DIM , SmoothingKernel kType >
OldKernel<DIM,kType>::OldKernel(	const int sle,
							const double hsml  ) :
	m_sle ( static_cast<SmoothingLengthEvolution> ( sle ) ),
	m_hsml( hsml ) {
}

// Std. Constructor
template < uint8_t DIM , SmoothingKernel kType >
OldKernel<DIM,kType>::OldKernel( 	const SmoothingLengthEvolution sle,
							const double hsml ) :
	m_sle ( sle ),
	m_hsml( hsml ) {

}

template < uint8_t DIM , SmoothingKernel kType >
OldKernel<DIM,kType>::~OldKernel() = default;

template < uint8_t DIM , SmoothingKernel kType >
template <typename T>
T OldKernel<DIM,kType>::get( const kMember & Member ) const{

	switch ( Member ){

		case kMember::DIM:
			return DIM;

		case kMember::kType:
			return kType;

		case kMember::sle:
			return m_sle;

		case kMember::hsml:
			return m_hsml;

		default:
			// Error

	}
}


// Specialized Kernels
// ---------------------------------------------------------------------------------------------

//template < uint8_t DIM >
//constexpr double OldKernel<DIM,SmoothingKernel::CubicSpline>::W() const {
//
//	const double f		= factor();
//	constexpr double q 	= .0;
//
//	return f * ( 2.0 / 3.0 - pow( q , 2.0 ) + 0.5 * pow( q , 3.0 ) );
//}

//template < uint8_t DIM >
//double OldKernel<DIM,SmoothingKernel::Gauss>::W() const {
//
//	const double f = 1.0 / ( pow( m_hsml, static_cast<double>( DIM ) ) * pow( CONST::PI , static_cast<double>( DIM ) / 2.0 ) );
//
//	return f;
//}

//template < uint8_t DIM >
//double OldKernel<DIM,SmoothingKernel::Quintic>::W() const {
//
//	const double f 		= factor();
//	constexpr double q	= .0;
//
//	return f * ( pow( 3.0 - q , 5.0 ) - 6.0 * pow( 2.0 - q , 5.0 ) + 15.0 * pow( 1.0 - q , 5.0 ) );
//
//}



//template < uint8_t DIM >
//double OldKernel<DIM,SmoothingKernel::CubicSpline>::W( const double r ) const {
//
//	const double f	= factor();
//	const double q 	= r / m_hsml;
//
//	double W;
//	if ( q >= 0.0 && q <= 1.0 ){
//		W = f * ( 2.0 / 3.0 - pow( q , 2.0 ) + 0.5 * pow( q , 3.0 ) );
//	}
//	else if ( q > 1.0 && q <= 2.0 ){
//		W = f * 1.0 / 6.0 * pow( 2.0 - q , 3.0 );
//
//	}
//	else {
//		W = .0;
//	}
//
//	return W;
//
//}

//template < uint8_t DIM >
//double OldKernel<DIM,SmoothingKernel::Gauss>::W( const double r ) const {
//
//	const double f = 1.0 / ( pow( m_hsml, static_cast<double>( DIM ) ) * pow( CONST::PI , static_cast<double>( DIM ) / 2.0 ) );
//	const double q = r / m_hsml;
//
//	double W;
//	if ( q >= .0 && q <= 3.0 ){
//		W = f * exp( -q * q );
//	}
//	else{
//		W = .0;
//	}
//
//	return W;
//
//}

//template < uint8_t DIM >
//double OldKernel<DIM,SmoothingKernel::Quintic>::W( const double r  ) const {
//
//	double W;
//	const double f 	= factor();
//	const double q	= r / m_hsml;
//
//	if( q >= .0 && q <= 1.0 ){
//		W = f * ( pow( 3.0 - q , 5.0 ) - 6.0 * pow( 2.0 - q , 5.0 ) + 15.0 * pow( 1.0 - q , 5.0 ) );
//	}
//	else if( q > 1.0 && q <= 2.0 ){
//		W = f * ( pow( 3.0 - q , 5.0 ) - 6.0 * pow( 2.0 - q , 5.0 ) );
//	}
//	else if( q > 2.0 && q <= 3.0 ){
//		W = f * pow( 3.0 - q , 5 );
//	}
//	else{
//		W = .0;
//	}
//
//	return W;
//}



//template < uint8_t DIM >
//constexpr v4df OldKernel<DIM,SmoothingKernel::CubicSpline>::dWdx() const {
//	return v4df( .0, .0, .0 , .0 );
//}

//template < uint8_t DIM >
//constexpr v4df OldKernel<DIM,SmoothingKernel::Gauss>::dWdx() const {
//	return v4df( .0, .0, .0 , .0 );
//}
//
//template < uint8_t DIM >
//constexpr v4df OldKernel<DIM,SmoothingKernel::Quintic>::dWdx() const {
//	return v4df( .0, .0, .0 , .0 );
//}
//


//template < uint8_t DIM >
//v4df OldKernel<DIM,SmoothingKernel::CubicSpline>::dWdx( const double r , const v4df & dx ) const {
//
//	const double f 	= factor();
//	const double q 	= r / m_hsml;
//
//
//	kOutput K;
//
//	v4df dWdx;
//	for( int d = 0; d < DIM; d++ ){
//		if ( q > 0.0 && q <= 1.0 ){
//			K.dwdx[d] = f * ( -2.0 + 3.0/2.0*q ) / pow( m_hsml , 2.0 ) * dx[d];
//		}
//		else if ( q > 1.0 && q <= 2.0 ){
//			K.dwdx[d] = - f * 1.0 / 6.0 * 3.0 * pow( 2.0 - q , 2.0 ) / m_hsml * ( dx[d] / r );
//		}
//		else {
//			K.dwdx[d] = .0;
//		}
//	}
//
//	return dWdx;
//
//
//}

//template < uint8_t DIM >
//v4df OldKernel<DIM,SmoothingKernel::Gauss>::dWdx( const double r , const v4df & dx ) const {
//
//	const double f 	= 1.0 / ( pow( m_hsml, static_cast<double>( DIM ) ) * pow( CONST::PI , static_cast<double>( DIM ) / 2.0 ) );
//	const double q	= r / m_hsml;
//
//	v4df dWdx;
//	const double w 	= W( r );
//	for( int d = 0; d < DIM; d++ ){
//		if ( q > .0 && q <= 3.0 ){
//			dWdx[d] = w * ( - 2.0 * dx[d] / pow( m_hsml , 2.0 ) );
//		}
//		else{
//			dWdx[d] = .0;
//		}
//	}
//
//	return dWdx;
//}

//template < uint8_t DIM >
//v4df OldKernel<DIM,SmoothingKernel::Quintic>::dWdx( const double r , const v4df & dx ) const {
//
//
//	v4df dWdx;
//
//	const double f 	= factor();
//	const double q	= r / m_hsml;
//
//
//	for( int d = 0; d < DIM; d++ ){
//		if( q > .0 && q <= 1.0 ){
//			dWdx[d] = f * ( ( -120.0 + 120.0*q - 50.0 * pow( q , 2.0 ) ) / pow( m_hsml , 2.0 ) * dx[d] );
//		}
//		else if( q > 1.0 && q <= 2.0 ){
//			dWdx[d] = f * ( -5.0 * pow( 3.0 - q , 4.0 ) + 30.0 * pow( 2.0 - q , 4.0 ) )/ m_hsml * ( dx[d] / r );
//		}
//		else if( q > 2.0 && q <= 3.0 ){
//			dWdx[d] = f * ( -5.0 * pow( 3.0 - q , 4.0 ) ) / m_hsml * ( dx[d] / r );
//		}
//		else{
//			dWdx[d] = .0;
//		}
//	}
//
//	return dWdx;
//}

////  Cubic Spline Factors
//double OldKernel<1,SmoothingKernel::CubicSpline>::factor () const {
//	return 1.0 / m_hsml;
//}
//
//double OldKernel<2,SmoothingKernel::CubicSpline>::factor () const {
//	return 15.0 /( 7.0 * CONST::PI * pow( m_hsml , 2.0 ) );
//}
//
//double OldKernel<3,SmoothingKernel::CubicSpline>::factor () const {
//	return 3.0 / ( 2.0 * CONST::PI * pow( m_hsml , 3.0 ) );
//}


//// Quintic Factors
// double OldKernel<1,SmoothingKernel::Quintic>::factor () const {
//	return 1.0 / ( 120.0 * m_hsml );
//}
//
//double OldKernel<2,SmoothingKernel::Quintic>::factor () const {
//	return 7.0 / ( 478.0 * CONST::PI * pow( m_hsml , 2.0 ) );
//}
//
//double OldKernel<3,SmoothingKernel::Quintic>::factor () const {
//	return factor = 1.0 / ( 120.0 * CONST::PI * pow( m_hsml , 3.0 ) );
//}


} /* namespace SPH */






























































