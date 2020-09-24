/*
 * SPH.h
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */




#ifndef SPH_SPHSIM_H_
#define SPH_SPHSIM_H_

#include <math.h>
#include <iostream>
#include <sstream>

#include "Primitives/DataTypes/MyTypes.h"
#include "Utilities/Timer.h"
#include "Settings/Settings.h"
#include "Particles/Particle.h"
#include "Particles/InteractionPair.h"
#include "Particles/NeighbourSearch.h"
#include "IO/Input.h"
#include "IO/Output.h"


namespace SPH {


template<uint8_t DIM>
class SimulationBox  {

private:

	const vDdf<DIM> m_Min;	//	xMin	: Lower limit of allowed x-regime	   -10.0
	const vDdf<DIM> m_Max;	// 	xMax	: Upper limit of allowed x-regime		10.0

public:

	SimulationBox( const Settings::Dictionary & s );

	virtual ~SimulationBox() = default;

};

template < uint8_t DIM >
SimulationBox<DIM>::SimulationBox( const Settings::Dictionary & s ) :
	m_Min( s.getAsVector<double>("origin") ),
	m_Max( m_Min + s.getAsVector<double>("extent") ){

}

template < uint8_t DIM >
class sphSim {

	typedef Particle<DIM,ParticleType::Physical> 	physPart;
	typedef Particle<DIM,ParticleType::Virtual> 	virtPart;

	typedef InteractionPair<DIM,ParticleType::Physical,ParticleType::Physical> 	ppInt;
	typedef InteractionPair<DIM,ParticleType::Physical,ParticleType::Virtual> 	pvInt;

	typedef NeighbourSearch<DIM,ParticleType::Physical>	SelfSearch_t;
	typedef std::unique_ptr<SelfSearch_t> 				SelfSearchPtr;

	typedef NeighbourSearch<DIM,ParticleType::Physical,ParticleType::Virtual>	pvSearch_t;
	typedef std::unique_ptr<pvSearch_t> 										pvSearchPtr;

private:

	const Settings::Settings 			m_Settings;

	int m_CurrentTS;

	double m_dt;											// 	dt 		: time step [s]
	double m_time; 											// 	time	: current simulation time [s]
	double m_tSim; 											// 	tSim 	: End time of simulation [s]



	std::vector<physPart> 		m_Particles;
	std::vector<virtPart> 		m_virtParticles;

	SelfSearchPtr 				m_pS;
	pvSearchPtr 				m_pvS;

	std::vector<ppInt> 			m_Interactions;
	std::vector<pvInt> 			m_virtInteractions;
	
	SimulationBox<DIM> 			m_SimulationBox;

	Output<DIM,ParticleType::Physical> m_outParticles;


public:

	// Empty Constructor
	sphSim();

	// Std. Constructor
	sphSim( const Settings::Settings & s);

	// Std. Destructor
	virtual ~sphSim();

	// run Simulation
	void run();

	// Time Integration
	void time_integration();

	// Calculate RHS
	void single_step(	vNDdf<DIM> 	& dvdt,
						vNdf 		& dudt,
						vNdf 		& dsdt,
						vNdf 		& tdsdt,
						vNdf		& drhodt,
						vNDdf<DIM>	& av );

	//template < DensityCalculation DC >
	void comp_density ( vNdf & rho );

	// Summation density
	void sum_density ( vNdf & rho );

	void cspm_density ( vNdf & rho );

	// Density with continuity approach
	void con_density ( vNdf & rho );

	// Calculate internal forces
	void int_force(	vNdf 		& tdsdt,
					vNdf 		& dedt,
					vNDdf<DIM> 	&	dvdt ) const;


	// Calculate the artificial viscosity
	void art_visc(	vNDdf<DIM> 	& dvdt,
					vNdf  		& dedt);

	// Calculate the artificial heat
	void art_heat( vNdf  & dedt );

	// Calculate the external forces
	void ext_force( vNDdf<DIM> & dvdt );

	// Evolve the smoothing length
	void h_upgrade( const vNdf & drhodt  );

	// 	Subroutine to calculate the average velocity to correct velocity
	// 	for preventing penetration (Monaghan, 1992)
	void av_vel( vNDdf<DIM> & av );


};

} /* namespace SPH */

#endif /* SPH_SPHSIM_H_ */


















