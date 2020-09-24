/*
 * Particle.h
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */




#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <math.h>
#include <cmath>
#include <math.h>
#include <array>

#include "Kernels/Kernels.h"
#include "Primitives/DataTypes/MyTypes.h"
#include "Settings/Settings.h"


namespace SPH {

enum class pMember {
	DIM,			//	DIM 	: Dimension
	x,				// 	x		: coordinates of particles
	v,				// 	v		: velocities of particles
	mass,			// 	mass	: mass of particles
	rho,			//	rho		: dnesities of particles
	p,				//	p		: pressure of particles
	u,				//	u		: internal energy of particles
	c,				//	c		: sound velocity of particles
	s,				//	s		: entropy of particles
	e,				//	e		: total energy of particles
	W,				//	W 		: Kernel function of Particle
	iType,			// 	itype	: types of particles
	pType			// 	pType	: physical, virtual Particle
};



enum class artEOS {
	NONE,							// 	0	: Not Specified
	Monaghan,						//	1	: Monaghan, 1994
	Morris							//	2	: Morris, 	1997
};

template < uint8_t DIM , ParticleType iType , ParticleType jType  >
class InteractionPair;

template<uint8_t DIM , ParticleType pType >
class Particle {

	typedef Kernel<DIM> 				Kernel_t;
	typedef std::unique_ptr<Kernel_t>	KernelPtr;
private:

	vDdf<DIM> m_x;						// 	x		: coordinates of particles
	vDdf<DIM> m_v;						// 	v		: velocities of particles
	vDdf<DIM> m_a; 						// 	a 		: dv/dt acceleration of particle

	double m_mass;					// 	mass	: mass of particles
	double m_rho;					//	rho		: dnesities of particles
	double m_mByRho;
	double m_p;						//	p		: pressure of particles
	double m_u;						//	u		: internal energy of particles
	double m_c;						//	c		: sound velocity of particles
	double m_s;						//	s		: entropy of particles
	double m_e;						//	e		: total energy of particles

	KernelPtr m_K;					//	W 		: Kernel function of Particle

	FluidType m_fType;				// 	itype	: types of particles


public:

	// Empty constructor
	Particle();

	// Copy constructor
	Particle( const Particle & p );

	// Normal constructor
	Particle(	const vDdf<DIM> & x,
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
				const FluidType itype );

	// Normal constructor
	Particle(	const vDdf<DIM> & x,
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
				const uint8_t itype );


	// Std. destructor
	virtual ~Particle();

	bool isPhysical() const;

	bool isVirtual() const;

	void set_acc( vDdf<DIM> a );

	// 	Subroutine to define the fluid particle viscosity
	const double viscosity() const;

	// 	Subroutine to calculate the smoothing kernel wij and its
	// 	derivatives dwdxij.
	void kernel(const double r,
				const vDdf<DIM> & dx,
				double & W,
				vDdf<DIM> & dWdx ) const;

	double hsml() const;

	double k() const;

	// Equation of state
	void eos();

	// 	Gamma law EOS: subroutine to calculate the pressure and sound
	void p_gas(	);

	// 	Artificial equation of state for the artificial compressibility
	void p_art_water();

	symdf<DIM> update_tMatrix( const symdf<DIM> & h ) const;

	double vcc ( const double hvcc ) const;

	double update_tdsdt( const symdf<DIM> & t  ) const;

	double p() const;

	double rho() const;

	double m ( ) const;

	double mByRho () const;

	double self_rho () const;

	double self_W () const;

	double update_u (	const double dt,
						const Symmetry sym,
						const double dudt );

	void update_u (		const double dt,
						const Symmetry sym,
						const double dudt,
						const double u_min );

	double update_rho(	const double dt,
						const Symmetry sym,
						const double drhodt );

	void update_rho(	const double dt,
						const Symmetry sym,
						const double drhodt,
						const double rho_min );

	void set_rho( 		const double rho );

	vDdf<DIM> update_v (const double dt,
						const vDdf<DIM> & av = vDdf<DIM>() );

	void update_v (		const double dt,
						const vDdf<DIM> & av,
						const vDdf<DIM> & v_min );


	void update_x (	const double dt );

	void update_h( const double drhodt , const double dt );

	void h_densityAlgebraic( const double fac );

	void h_smlODE( const double drhodt , const double dt );

	template<uint8_t D , ParticleType P >
	friend std::ostream & operator<< ( std::ostream & s, const Particle<D,P> & p );

	template<uint8_t D , ParticleType P1, ParticleType P2 >
	friend vDdf<D> dx(const Particle<D,P1> & p1, const Particle<D,P2> & p2 );

	template<uint8_t D , ParticleType P1, ParticleType P2 >
	friend vDdf<D> dv(const Particle<D,P1> & p1, const Particle<D,P2> & p2);

	template <typename T>
	T get( const pMember & Member ) const;

	template < uint8_t D , ParticleType I , ParticleType J  >
	friend class InteractionPair;

	template < uint8_t D, ParticleType Type >
	friend class Output;
};

template<uint8_t D , ParticleType P1, ParticleType P2 >
inline
vDdf<D> dx(const Particle<D,P1> & p1, const Particle<D,P2> & p2 ) {
	return p1.m_x - p2.m_x;
}

template<uint8_t D , ParticleType P1, ParticleType P2 >
inline
vDdf<D> dv(const Particle<D,P1> & p1, const Particle<D,P2> & p2) {
	return p1.m_v - p2.m_v;
}


//template< ParticleType pType >
//class Particle<3,pType> { double tdsdt ( const symdf<3> & t ) const; };

} /* namespace SPH */

#endif /* PARTICLE_H_ */




























