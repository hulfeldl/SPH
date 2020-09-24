/*
 * InteractionPair.h
 *
 *  Created on: Aug 10, 2020
 *      Author: hulfeldl
 */




#ifndef PARTICLES_INTERACTIONPAIR_H_
#define PARTICLES_INTERACTIONPAIR_H_

#include <array>

#include "Primitives/DataTypes/MyTypes.h"
#include "Particles/Particle.h"

namespace SPH {

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
class InteractionPair {

private:

	const int m_i;
	const int m_j;

	const Particle<DIM,iType> & m_iParticle;
	const Particle<DIM,jType> & m_jParticle;

	double 		m_W;
	vDdf<DIM> 	m_dWdx;

	double 		m_r;
	vDdf<DIM> 	m_dx;
	vDdf<DIM> 	m_dv;


public:

	InteractionPair();

	InteractionPair( const int i, const int j );

	InteractionPair(const int i,
					const int j,
					const double r,
					const vDdf<DIM> & dx,
					const Particle<DIM,iType> & iParticle,
					const Particle<DIM,jType> & jParticle );

	virtual ~InteractionPair();

	int i() const;
	int j() const;

	double 	W() 	const;
	vDdf<DIM> 	dWdx() 	const;
	vDdf<DIM> 	dv() 	const;

	symdf<DIM> hMatrix() const;

	double hvcc() const;

	void calc_wi( vNdf & wi );

	void calc_rho( vNdf & rho );

	void calc_con( vNdf & drhodt );


	void calc_tMatrix(	vNSdf<DIM> 	& tMatrix,
						vNdf  		& vcc ) const;

	void SPH_Algorithm1 ( 	const vNSdf<DIM> 	& t,
							const bool isViscous,
							vNdf  				& dedt,
							vNDdf<DIM> 			& dvdt ) const;

	void SPH_Algorithm2 ( 	const vNSdf<DIM> 	& t,
							const bool isViscous,
							vNdf  				& dedt,
							vNDdf<DIM>			& dvdt ) const;

	void art_viscosity(	vNDdf<DIM> & dvdt,
						vNdf  & dedt );

	void calc_vcc ( vNdf & vcc );

	void art_heat ( 	const vNdf & vcc,
											vNdf & dedt );

	void boundary_force( vNDdf<DIM> & dvdt );

	vDdf<DIM> int_phys_virt();

	void av_vel ( vNDdf<DIM> & av ) const;

};

} /* namespace SPH */

#endif /* PARTICLES_INTERACTIONPAIR_H_ */



















