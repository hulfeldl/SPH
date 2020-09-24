/*
 * InteractionPair.cpp
 *
 *  Created on: Aug 10, 2020
 *      Author: hulfeldl
 */

#include "InteractionPair.h"

namespace SPH {

//template < uint8_t DIM , ParticleType iType , ParticleType jType  >
//InteractionPair<DIM,iType,jType>::InteractionPair() {
//	// TODO Auto-generated constructor stub
//
//}
//
//template < uint8_t DIM , ParticleType iType , ParticleType jType  >
//InteractionPair<DIM,iType,jType>::InteractionPair( const int i, const int j ){
//
//}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
InteractionPair<DIM,iType,jType>::InteractionPair(	const int i,
													const int j,
													const double r,
													const vDdf<DIM> & dx,
													const Particle<DIM,iType> & iParticle,
													const Particle<DIM,jType> & jParticle ) :
	m_i ( i ),
	m_j ( j ),
	m_iParticle ( iParticle ),
	m_jParticle ( jParticle ),
	m_r ( r ),
	m_dx ( dx ),
	m_dv ( jParticle.m_v - iParticle.m_v ) {

	double iW; vDdf<DIM> idWdx;
	iParticle.kernel( m_r, m_dx, iW, idWdx );

	double jW; vDdf<DIM> jdWdx;
	jParticle.kernel( m_r, m_dx, jW, jdWdx );

	m_W 	= 0.5 * ( iW + jW );
	m_dWdx 	= 0.5 * ( idWdx + jdWdx );
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
InteractionPair<DIM,iType,jType>::~InteractionPair() = default;

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
int InteractionPair<DIM,iType,jType>::i() const{
	return m_i;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
int InteractionPair<DIM,iType,jType>::j() const{
	return m_j;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
double 	InteractionPair<DIM,iType,jType>::W() const{
	return m_W;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
vDdf<DIM>	InteractionPair<DIM,iType,jType>::dWdx() const{
	return m_dWdx;

}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
vDdf<DIM>	InteractionPair<DIM,iType,jType>::dv() 	const {
	return m_dv;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
symdf<DIM> InteractionPair<DIM,iType,jType>::hMatrix() const {

	return h( m_dv, m_dWdx );
}


template < uint8_t DIM , ParticleType iType , ParticleType jType  >
double InteractionPair<DIM,iType,jType>::hvcc() const {

	double hvcc = m_dv ^ m_dWdx;
	return hvcc;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::calc_wi( vNdf & wi ) {

	wi[m_i] += 	m_jParticle.mByRho() * m_W;
	wi[m_j] += 	m_iParticle.mByRho() * m_W;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::calc_rho( vNdf & rho ) {

	rho[m_i] += m_jParticle.m() * m_W;
	rho[m_j] += m_iParticle.m() * m_W;
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::calc_con( vNdf & drhodt ) {

	double vcc = - ( m_dv ^ m_dWdx );

	drhodt[m_i] += m_jParticle.m() * vcc;
	drhodt[m_j] += m_iParticle.m() * vcc;

}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::calc_tMatrix(vNSdf<DIM> 	& tMatrix,
													vNdf  		& vcc ) const {

	const symdf<DIM> h 	= hMatrix();

	tMatrix[m_i] += m_jParticle.update_tMatrix( h );
	tMatrix[m_j] += m_iParticle.update_tMatrix( h );

	// Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz:
	double HVcc =  m_dv ^ m_dWdx;

	vcc[m_i] += m_jParticle.mByRho() * HVcc;
	vcc[m_j] += m_iParticle.mByRho() * HVcc;
}


template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::SPH_Algorithm1 (	const vNSdf<DIM> 	& t,
														const bool isViscous,
														vNdf  				& dedt,
														vNDdf<DIM>			& dvdt ) const {

	const double etai			= m_iParticle.viscosity();
	const double etaj 			= m_jParticle.viscosity();

	const double rhoij = 1.0 / ( m_iParticle.rho() * m_jParticle.rho() );

	vDdf<DIM> h = -( m_iParticle.p() + m_jParticle.p() ) * m_dWdx;
	if (isViscous) {
		h += ( etai * t[m_i] + etaj * t[m_j] ) * m_dWdx;
	}

	h 	*= rhoij;
	dvdt[m_i] += m_jParticle.m() * h;
	dvdt[m_j] -= m_iParticle.m() * h;

	double he = h ^ m_dv;

	dedt[m_i] += m_jParticle.m() * he;
	dedt[m_j] += m_iParticle.m() * he;


}


template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::SPH_Algorithm2 ( const vNSdf<DIM> 	& t,
														const bool isViscous,
														vNdf  				& dedt,
														vNDdf<DIM> 				& dvdt ) const {

	const double etai			= m_iParticle.viscosity();
	const double etaj 			= m_jParticle.viscosity();

	const double rhoi2 			= m_iParticle.rho() * m_iParticle.rho();
	const double rhoj2 			= m_jParticle.rho() * m_jParticle.rho();

	vDdf<DIM> h = -( m_iParticle.p() / rhoi2 + m_jParticle.p() / rhoj2 ) * m_dWdx;
	if (isViscous) {
		h += ( (etai/ rhoi2)  * t[m_i] + (etaj / rhoj2) * t[m_j] ) * m_dWdx;
	}

//	std::cout << "::SPH_Algorithm2 " << "dWdx: " << m_dWdx << std::endl;
//	std::cout << "::SPH_Algorithm2 " << "h: " << h << std::endl;
//	std::cout << "::SPH_Algorithm2 " << "dvdt_i: " << dvdt[m_i] << "dvdt_j: " << dvdt[m_j] << std::endl;
//
	dvdt[m_i] += m_jParticle.m() * h;
	dvdt[m_j] -= m_iParticle.m() * h;

//	std::cout << "::SPH_Algorithm2 " << "dvdt_i: " << dvdt[m_i] << "dvdt_j: " << dvdt[m_j] << std::endl;

	double he = h ^ m_dv;
	dedt[m_i] += m_jParticle.m() * he;
	dedt[m_j] += m_iParticle.m() * he;

}


template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::art_viscosity(	vNDdf<DIM> 	& dvdt,
														vNdf  		& dedt ) {

	// Parameter for the artificial viscosity:
	// Shear viscosity
	constexpr double alpha 	= 1.0;

	// Bulk viscosity
	constexpr double beta 	= 1.0;

	// Parameter to avoid singularities
	constexpr double etq	= 0.1;

	const double mhsml = 0.5 * ( m_iParticle.hsml() + m_jParticle.hsml() );

	double vr = - ( m_dv ^ m_dx );
	double rr = m_dx ^ m_dx;

	// Artificial viscous force only if v_ij * r_ij < 0
	if (vr < .0){

		// Calculate muv_ij = hsml v_ij * r_ij / ( r_ij A 2 + hsml*2 etq*2 )
		const double muv = mhsml * vr / ( rr + mhsml * mhsml * etq * etq );

		// Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij*2) / rho_ij
		const double mc 	= 0.5 * ( m_iParticle.m_c + m_jParticle.m_c );
		const double mrho 	= 0.5 * ( m_iParticle.m_rho + m_jParticle.m_rho );
		const double piv	= ( beta * muv - alpha * mc ) * muv / mrho;


		// Calculate SPH sum for artificial viscous force
		const vDdf<DIM> h	= -piv * m_dWdx;
		dvdt[m_i]   += m_jParticle.m_mass * h;
		dvdt[m_j]	-= m_iParticle.m_mass * h;

		const double he = h ^ m_dv;
		dedt[m_i] 	+= m_jParticle.m_mass * he;
		dedt[m_j] 	+= m_iParticle.m_mass * he;
	}
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::calc_vcc ( vNdf & vcc ) {

	const double HVCC = hvcc();

	vcc[m_i] += m_jParticle.m_mByRho * HVCC;
	vcc[m_j] += m_iParticle.m_mByRho * HVCC;
}


template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::art_heat ( 	const vNdf & vcc,
													vNdf & dedt ) {

	// Parameter for the artificial heat conduction:
	constexpr double g1	= 0.1;
	constexpr double g2	= 1.0;

	const double mhsml 	= 0.5 * m_iParticle.hsml() + m_jParticle.hsml();
	const double mrho	= 0.5 * ( m_iParticle.m_rho + m_jParticle.m_rho );

	double rdWdx 	= m_dx ^ m_dWdx;
	double rr 		= m_dx ^ m_dx;

	const double ihsml 	= m_iParticle.hsml();
	const double jhsml 	= m_jParticle.hsml();
	const double mui	= g1 * ihsml * m_iParticle.m_c + g2 * pow( ihsml , 2.0 ) * ( fabs( vcc[m_i]) - vcc[m_i] );
	const double muj	= g1 * jhsml * m_jParticle.m_c + g2 * pow( jhsml , 2.0 ) * ( fabs( vcc[m_j]) - vcc[m_j] );
	const double muij	= 0.5 * ( mui + muj );

	const double h 		= muij / ( mrho * ( rr + 0.01 * mhsml*mhsml ) ) * rdWdx;

	const double iu 	= m_iParticle.m_u;
	const double ju 	= m_jParticle.m_u;

	dedt[m_i] += m_jParticle.m_mass * h * ( iu - ju );
	dedt[m_j] += m_iParticle.m_mass * h * ( ju - iu );
}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::boundary_force(	vNDdf<DIM> & dvdt ) {


	if ( iType == ParticleType::Physical && jType == ParticleType::Virtual ){
		dvdt[m_i] += int_phys_virt();
	}
	else if ( iType == ParticleType::Virtual && jType == ParticleType::Physical ){
		dvdt[m_j] -= int_phys_virt();
	}




}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
vDdf<DIM> InteractionPair<DIM,iType,jType>::int_phys_virt() {

	// Boundary particle force and penalty anti-penetration force.
	constexpr double rrO 	= 1.25e-5;
	constexpr double dd	 	= 1.0e-2;
	constexpr double p1 	= 12.0;
	constexpr double p2 	= 4.0;

	vDdf<DIM> dvdt{};
	if( m_r < rrO ){
		const double f = ( pow( rrO / m_r , p1 ) - pow( rrO / m_r , p2 ) ) / ( pow( m_r , 2.0 ) );

		for( int d = 0; d < DIM; d++ ){
			dvdt[d] = dd * m_dx[d] * f;
		}
	}

	return dvdt;

}

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
void InteractionPair<DIM,iType,jType>::av_vel ( vNDdf<DIM> & av ) const {

	const double mrho 	= 0.5 * ( m_iParticle.rho() + m_jParticle.rho() );

	av[m_i] += m_jParticle.m() * ( m_W / mrho ) * m_dv;
	av[m_j] -= m_iParticle.m() * ( m_W / mrho ) * m_dv;
}

template class InteractionPair<1,ParticleType::Physical,ParticleType::Physical>;
template class InteractionPair<2,ParticleType::Physical,ParticleType::Physical>;
template class InteractionPair<3,ParticleType::Physical,ParticleType::Physical>;

template class InteractionPair<1,ParticleType::Physical,ParticleType::Virtual>;
template class InteractionPair<2,ParticleType::Physical,ParticleType::Virtual>;
template class InteractionPair<3,ParticleType::Physical,ParticleType::Virtual>;

} /* namespace SPH */


























