/*
 * Particle.cpp
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */




#include <sstream>

#include "Particle.h"


namespace SPH {

template<uint8_t DIM , ParticleType pType >
Particle<DIM,pType>::Particle() :
	m_x (),
	m_v (),
	m_a (),
	m_mass ( .0 ),
	m_rho ( .0 ),
	m_mByRho( m_mass/m_rho ),
	m_p ( .0 ),
	m_u ( .0 ),
	m_c ( .0 ),
	m_s ( .0 ),
	m_e ( .0 ),
	m_K (),
	m_fType ( FluidType::NONE ) {

}

// Copy constructor
template<uint8_t DIM , ParticleType pType >
Particle<DIM,pType>::Particle( const Particle & p ) :
	m_x ( p.m_x ),
	m_v ( p.m_v ),
	m_a (),
	m_mass ( p.m_mass ),
	m_rho ( p.m_rho ),
	m_mByRho( m_mass/m_rho ),
	m_p ( p.m_p ),
	m_u ( p.m_u ),
	m_c ( p.m_c ),
	m_s ( p.m_s ),
	m_e ( p.m_e ),
	m_K ( p.m_K->clone() ),
	m_fType ( p.m_fType ) {

}

// Normal constructor
template<uint8_t DIM , ParticleType pType >
Particle<DIM,pType>::Particle(	const vDdf<DIM> & x,
								const vDdf<DIM> & v,
								const double m,
								const double rho,
								const double p,
								const double u,
								const double c,
								const double s,
								const double e,
								const SmoothingKernel kType,
								const double hsml,
								const FluidType itype) :
	m_x ( x ),
	m_v ( v ),
	m_a (),
	m_mass ( m ),
	m_rho ( rho ),
	m_mByRho( m_mass/m_rho ),
	m_p ( p ),
	m_u ( u ),
	m_c ( c ),
	m_s ( s ),
	m_e ( e ),
	m_K ( std::make_unique<CubicSplineKernel<DIM>>(SmoothingLengthEvolution::DensityAlgebraic , hsml) ),
	m_fType ( itype ) {

}

// Normal constructor
template<uint8_t DIM , ParticleType pType >
Particle<DIM,pType>::Particle(	const vDdf<DIM> & x,
								const vDdf<DIM> & v,
								const double m,
								const double rho,
								const double p,
								const double u,
								const double c,
								const double s,
								const double e,
								const uint8_t kType,
								const double hsml,
								const uint8_t itype) :
	m_x ( x ),
	m_v ( v ),
	m_a (),
	m_mass ( m ),
	m_rho ( rho ),
	m_mByRho( m_mass/m_rho ),
	m_p ( p ),
	m_u ( u ),
	m_c ( c ),
	m_s ( s ),
	m_e ( e ),
	m_K ( std::make_unique<CubicSplineKernel<DIM>>(SmoothingLengthEvolution::DensityAlgebraic , hsml) ),
	m_fType ( static_cast<FluidType>(itype) ) {

}

template<uint8_t DIM , ParticleType pType >
Particle<DIM,pType>::~Particle() = default;

template<uint8_t DIM , ParticleType pType >
bool Particle<DIM,pType>::isPhysical() const {
	return pType == ParticleType::Physical;
}

template<uint8_t DIM , ParticleType pType >
bool Particle<DIM,pType>::isVirtual() const {
	return pType == ParticleType::Virtual;
}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::set_acc( vDdf<DIM> a ) {
	m_a = a;
}

//----------------------------------------------------------------------------
// 	Subroutine to define the fluid particle viscosity

//	itype	: Type of particle				[in]
//	eta		: Dynamic viscosity				[out]
//----------------------------------------------------------------------------
template<uint8_t DIM , ParticleType pType >
const double Particle<DIM,pType>::viscosity() const {

	double eta;

	switch( m_fType ){

		case FluidType::IdealGas:
			eta = .0;
			break;

		case FluidType::Water:
			eta = 1.0e-3;
			break;

		case FluidType::NONE:
			MYASSERT (false, "Particle has no Fluid Type! Unable to update Equation of State!");
	}

	return eta;

}

//----------------------------------------------------------------------------
// 	Subroutine to calculate the smoothing kernel wij and its
// 	derivatives dwdxij.

// 	r 		: Distance between particles i and j												[in]
// 	dx		: x-, y- and z-distance between i and j 											[in]
//	kOutput	: Kernel & Derivative of kernel with respect to x, y and z for interaction pair 	[out]
//----------------------------------------------------------------------------
template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::kernel(	const double r,
									const vDdf<DIM> & dx,
									double & W,
									vDdf<DIM> & dWdx ) const {
	W 		= m_K->W( r );
	dWdx	= m_K->dWdx( r , dx );
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::hsml() const {
	return m_K->hsml();
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::k() const {
	return m_K->k();
}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::eos(){

	switch ( m_fType ) {

		case FluidType::IdealGas:
			p_gas();
			break;

		case FluidType::Water:
			p_art_water();
			break;

		case FluidType::NONE:
			MYASSERT (false, "Particle has no Fluid Type! Unable to update Equation of State!");

	}

}


//----------------------------------------------------------------------------
// 	Gamma law EOS: subroutine to calculate the pressure and sound

//	rho	: Density 			[in]
//	u	: Internal energy	[in]
//	p	: Pressure			[out]
//	c	: sound velocity	[out]
//----------------------------------------------------------------------------
template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::p_gas(){

	// For air (ideal gas)
	constexpr double gamma	= 1.4;
	m_p					= (gamma - 1)*m_rho * m_u;
	m_c 				= sqrt( ( gamma - 1 ) * m_u );

}

//----------------------------------------------------------------------------
// 	Artificial equation of state for the artificial compressibility

// 	rho : Density 			[in]
// 	u 	: Internal energy 	[in]
// 	p 	: Pressure 			[out]
// 	c 	: sound velocity 	[out]

// 	Equation of state for artificial compressibility
//----------------------------------------------------------------------------
template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::p_art_water(	){


	const artEOS eosType 	= artEOS::Monaghan;

	switch( eosType ){

		case artEOS::Monaghan: {

			// Artificial EOS, Form 1 (Monaghan, 1994)
			const double gamma	= 7.0;
			const double rho0	= 1000.0;
			const double b 		= 1.013e5;

			m_p 		= b*( pow( m_rho / rho0 , gamma ) - 1);
			m_c 		= 1480;
			break;
		}

		case artEOS::Morris:
			// Artificial EOS, Form 2 (Morris, 1997)
			m_c  		= 0.01;
			m_p 		= m_c * m_c * m_rho;
			break;

		case artEOS::NONE:
			MYASSERT (false,"No artificial Equation of State model chosen! \
							 Unable to calculate pressure for Water!");

	}

}

template<uint8_t DIM , ParticleType pType >
symdf<DIM> Particle<DIM,pType>::update_tMatrix( const symdf<DIM> & h ) const {
	return m_mByRho * h;
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::update_tdsdt( const symdf<DIM> & t ) const {
	return 0.5 * viscosity() / m_rho * t ^ t;
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::p() const {
	return m_p;
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::rho() const {
	return m_rho;
}


template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::m () const {
	return m_mass;
}


template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::mByRho () const {
	return ( m_mass / m_rho );
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::self_rho () const {
	return m_K->W() * m_mass;
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::self_W () const {
	return m_K->W() * mByRho();
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::update_u (	const double dt,
							const Symmetry sym,
							const double dudt ) {

	double u_min 	= m_u;
	double temp_u	= .0;

	// TODO: Wierd Symmetry thing, implement correct symmetry conditions
	if ( DIM == 1 ){
		temp_u	= -static_cast<int>( sym ) * m_p * m_v[0] / ( m_x[0] * m_rho );
	}

	m_u += dt * ( dudt + temp_u );

	if( m_u < 0 ){
		m_u = .0;
	}

	return u_min;
}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::update_u (	const double dt,
								const Symmetry sym,
								const double dudt,
								const double u_min ) {

	double temp_u	= .0;

	// TODO: Wierd Symmetry thing, implement correct symmetry conditions
	if ( DIM == 1 ){
		temp_u	= -static_cast<int>( sym ) * m_p * m_v[0] / ( m_x[0] * m_rho );
	}

	m_u = u_min + dt * ( dudt + temp_u );

	if( m_u < 0 ){
		m_u = .0;
	}
}

template<uint8_t DIM , ParticleType pType >
double Particle<DIM,pType>::update_rho(	const double dt,
									const Symmetry sym,
									const double drhodt ){

	double rho_min	= m_rho;
	double temp_rho	= .0;

	// TODO: Wierd Symmetry thing, implement correct symmetry conditions
	if ( DIM == 1 ){
		temp_rho = -static_cast<int>( sym ) * m_rho * m_v[0] / m_x[0];
	}

	m_rho += dt * ( drhodt + temp_rho );

	return rho_min;
}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::update_rho(	const double dt,
							const Symmetry sym,
							const double drhodt,
							const double rho_min ) {

	double temp_rho	= .0;

	// TODO: Wierd Symmetry thing, implement correct symmetry conditions
	if ( DIM == 1 ){
		temp_rho = -static_cast<int>( sym ) * m_rho * m_v[0] / m_x[0];
	}

	m_rho = rho_min + dt * ( drhodt + temp_rho );

}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::set_rho( 		const double rho ) {
	m_rho = rho;
}

template<uint8_t DIM , ParticleType pType >
vDdf<DIM> Particle<DIM,pType>::update_v (const double dt,
									const vDdf<DIM> & av ) {

	vDdf<DIM> v_min;

	for( int d = 0; d < DIM; d++ ){
		v_min[d]	= m_v[d];
		m_v[d] 	   += dt * m_a[d] + av[d];
	}

	return v_min;

}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::update_v (const double dt,
									const vDdf<DIM> & av,
									const vDdf<DIM> & v_min ) {

	for( int d = 0; d < DIM; d++ ){
		m_v[d] 		= v_min[d] + dt * m_a[d] + av[d];
	}

}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::update_x (	const double dt ) {

	for( int d = 0; d < DIM; d++ ){
		m_x[d] += dt * m_v[d];
	}
}

template<uint8_t DIM , ParticleType pType >
void Particle<DIM,pType>::update_h( const double drhodt , const double dt ) {

	switch ( m_K->sle() ){

	case SmoothingLengthEvolution::DensityAlgebraic: {
		const double hsml = 2.0 * ( pow( m_mass / m_rho , 1.0/static_cast<double>(DIM) ) );
		m_K->update_hsml( hsml );
		m_K->update_f();
		break;
	}

	case SmoothingLengthEvolution::SmoothingLengthODE: {
		// dh/dt = (-1/dim)*(h/rho)*(drho/dt).
		const double dhdt	= (-1.0 / static_cast<double>(DIM) ) * ( m_K->hsml() / m_rho ) * drhodt;
		m_K->update_hsml( dhdt , dt );
		m_K->update_f();
		break;
	}

	case SmoothingLengthEvolution::OtherApproach:

	case SmoothingLengthEvolution::NONE:
		break;
	}

}

//template<uint8_t DIM , ParticleType pType >
//template <typename T>
//T Particle<DIM,pType>::get( const pMember & Member ) const {
//
//
//	switch ( Member ){
//
//		case pMember::DIM:
//			return DIM;
//
//		case pMember::x:
//			return m_x;
//
//		case pMember::v:
//			return m_v;
//
//		case pMember::mass:
//			return m_mass;
//
//		case pMember::rho:
//			return m_rho;
//
//		case pMember::p:
//			return m_p;
//
//		case pMember::u:
//			return m_u;
//
//		case pMember::c:
//			return m_c;
//
//		case pMember::s:
//			return m_s;
//
//		case pMember::e:
//			return m_e;
//
//		case pMember::W:
//			return m_W;
//
//		case pMember::iType:
//			return m_fType;
//
//		case pMember::pType:
//			return m_pType;
//
//		default:
//			// Error
//	}
//
//}

template<uint8_t DIM , ParticleType pType >
std::ostream & operator<< ( std::ostream & os, const Particle<DIM,pType> & p ){

	os << "x:            " << p.m_x << std::endl;
	os << "v:            " << p.m_v << std::endl;
	os << "a:            " << p.m_a << std::endl;

	os << "rho:          " << p.m_rho << std::endl;
	os << "m/rho:        " << p.m_mByRho << std::endl;
	os << "p:            " << p.m_p << std::endl;
	os << "u:            " << p.m_u << std::endl;

	os << "hsml:         " << p.hsml() << std::endl;

	return os;

}

//template<uint8_t DIM , ParticleType pType >
//vDdf<DIM> Particle<DIM,pType>::dx(const Particle & jParticle) const {
//
//	v4df dx;
//	for ( int d = 0; d < DIM; d++ ) {
//		dx[d] 	= m_x - jParticle.m_x;
//	}
//
//	return dx;
//}
//
//template<uint8_t DIM , ParticleType pType >
//vDdf<DIM> Particle<DIM,pType>::dv(const Particle & jParticle) const {
//
//	v4df dv;
//	for( int d = 0; d < DIM; d++){
//		dv[d] 	= jParticle.m_v[d] - m_v[d];
//	}
//
//	return dv;
//}

template class Particle<1,ParticleType::Physical>;
template class Particle<2,ParticleType::Physical>;
template class Particle<3,ParticleType::Physical>;

template class Particle<1,ParticleType::Virtual>;
template class Particle<2,ParticleType::Virtual>;
template class Particle<3,ParticleType::Virtual>;

template std::ostream & operator<< ( std::ostream & os, const Particle<1,ParticleType::Physical> & p );
template std::ostream & operator<< ( std::ostream & os, const Particle<2,ParticleType::Physical> & p );
template std::ostream & operator<< ( std::ostream & os, const Particle<3,ParticleType::Physical> & p );

template std::ostream & operator<< ( std::ostream & os, const Particle<1,ParticleType::Virtual> & p );
template std::ostream & operator<< ( std::ostream & os, const Particle<2,ParticleType::Virtual> & p );
template std::ostream & operator<< ( std::ostream & os, const Particle<3,ParticleType::Virtual> & p );

} /* namespace SPH */





















