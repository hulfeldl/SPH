/*
 * SPH.cpp
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */

#include "sphSim.h"

namespace SPH {

template < uint8_t DIM >
sphSim<DIM>::sphSim( const Settings::Settings & s ) :
	m_Settings ( s ),
	m_CurrentTS (0),
	m_dt(m_Settings.dict().get<Settings::Dictionary>("time step control").get<double>("dt")),
	m_time(.0),
	m_tSim(m_dt*m_Settings.dict().get<Settings::Dictionary>("simulation").get<double>("n timesteps")),
	m_Particles (),
	m_virtParticles(),
	m_pS( 	std::make_unique<DirectSearch<DIM,ParticleType::Physical>>(m_Particles,m_Interactions,m_Settings.maxInteractions())),
	m_pvS( 	std::make_unique<DirectSearch<DIM,ParticleType::Physical,ParticleType::Virtual>>(m_Particles,m_virtParticles,m_virtInteractions,m_Settings.maxInteractions())),
	m_Interactions (),
	m_virtInteractions (),
	m_SimulationBox ( m_Settings.dict().get<Settings::Dictionary>("simulation box") ),
	m_outParticles ( m_Particles, m_Settings.opts().oDir(), m_Settings.outS().Fileformat()  ) {

	Input<DIM,ParticleType::Physical> PhysInput( m_Settings.opts().iDir()  );
	m_Particles = PhysInput.get_input( m_Settings.iniPhys(), m_Settings.pType()  );

	Input<DIM,ParticleType::Virtual> VirtInput( m_Settings.opts().iDir()  );
	m_virtParticles = VirtInput.get_input( m_Settings.iniVirt(), m_Settings.pType() );
}

template < uint8_t DIM >
sphSim<DIM>::~sphSim() = default;


// run Simulation
//----------------------------------------------------------------------------
// 	This is a three dimensional SPH code, the followings are the
// 	basic parameters needed in this code or calculated by this code

// 	mass	: mass of particles							[in]
//	ntotal	: total particle number						[in]
//	dt		: Time step used in the time integration	[in]
// 	itype	: types of particles						[in]
// 	x		: coordinates of particles					[in/out]
// 	vx		: velocities of particles					[in/out]
//	rho		: dnesities of particles					[in/out]
//	p		: pressure of particles						[in/out]
//	u		: internal energy of particles				[in/out]
//	hsml 	: smoothing lengths of particles			[in/out]
//	c		: sound velocity of particles				[out]
//	s		: entropy of particles						[out]
//	e		: total energy of particles					[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::run(){

	Timer timer;

	std::time_t s1, s2;

	timer.time_print();
	s1 = timer.time_elapsed();

//	bool runSimulation  = true;
//	while( runSimulation ){

		std::cout << "***************************************************" << std::endl;
		std::cout << "        Please input the maximal time steps 		 " << std::endl;
		std::cout << "***************************************************" << std::endl;


		time_integration( );

		std::string Filename = "Particles";
		m_outParticles.writeOutput( Filename, m_time );

//		std::cout << "***************************************************" << std::endl;
//		std::cout << "Are you going to run more time steps? (0=No, 1=yes)" << std::endl;
//		std::cout << "***************************************************" << std::endl;

//		std::cin >> runSimulation;
//	}

	timer.time_print();
	s2 = timer.time_elapsed();

	std::cout << "Elapsed CPU time = " << s2 - s1 << std::endl;

}


// Time Integration

//----------------------------------------------------------------------------

// 	x 			: coordinates of particles				[in/out]
// 	vx 			: velocities of particles				[in/out]
// 	mass 		: mass of particles						[in]
//	rho			: dnesities of particles				[in/out]
//	p			: pressure of particles					[in/out]
//	u			: internal energy of particles			[in/out]
//	c			: sound velocity of particles			[out]
//	s			: entropy of particles, not used here	[out]
//	e			: total energy of particles				[out]
//	itype		: types of particles					[in]
//					= 1 ideal gas
//					= 2 water
// 					= 3 tnt
// 	hsml		: smoothing lengths of particles		[in/out]
// 	ntotal		: total particle number					[in]
// 	maxtimestep	: maximum timesteps						[in]
//	dt			: timestep								[in]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::time_integration(){


	const size_t nParticles 	= m_Particles.size();

	vNDdf<DIM> v_min 	( nParticles, vDdf<DIM>() );
	vNDdf<DIM> dxdt 	( nParticles, vDdf<DIM>() );
	vNDdf<DIM> dvdt 	( nParticles, vDdf<DIM>() );
	vNDdf<DIM> av 		( nParticles, vDdf<DIM>() );

	vNdf u_min 		( nParticles, .0 );
	vNdf rho_min	( nParticles, .0 );
	vNdf dudt	 	( nParticles, .0 );
	vNdf drhodt 	( nParticles, .0 );
	vNdf dsdt 		( nParticles, .0 );
	vNdf t 			( nParticles, .0 );
	vNdf tdsdt 		( nParticles, .0 );

	const int save_step		= m_Settings.outS().SaveStep();
	Symmetry sym 			= Symmetry::NONE;

	std::cout << "**************************************************" << std::endl;
	std::cout << "Current number of time step = " << m_CurrentTS << "  	current time = " << m_time << std::endl;
	std::cout << "**************************************************" << std::endl;

	std::cout << std::endl;
	vNsu moniParticles = m_Settings.outS().MoniParticles();
	for ( size_t i = 0; i < moniParticles.size(); i++ ){
		const size_t iMoni 	= moniParticles[i];

		std::cout << "**** Information for particle **** " << iMoni << std::endl;
		std::cout << m_Particles[iMoni];
		std::cout << std::endl;

	}

	while ( m_time < m_tSim ) {

		// If not first time step, then update thermal energy, density and
		// velocity half a time step
		if ( m_CurrentTS != 0 ){

			for( size_t i = 0; i < nParticles; i++ ){

				u_min[i]	= m_Particles[i].update_u( 0.5 * m_dt, sym, dudt[i] );


				if ( m_Settings.get<DensityCalculation>() == DensityCalculation::ContinuityEquation ){
					rho_min[i] 	= m_Particles[i].update_rho( 0.5 * m_dt, sym, drhodt[i] );
				}

				v_min[i] = m_Particles[i].update_v( 0.5 * m_dt );
			}
		}

		// Definition of variables out of the function vector:
		single_step( dvdt, dudt, dsdt, tdsdt, drhodt, av );

		if ( m_CurrentTS == 0 ){

			for( size_t i = 0; i < nParticles; i++ ){

				m_Particles[i].update_u( 0.5 * m_dt, sym, dudt[i] );


				if ( m_Settings.get<DensityCalculation>() == DensityCalculation::ContinuityEquation ){
					m_Particles[i].update_rho( 0.5 * m_dt, sym, drhodt[i] );
				}

				m_Particles[i].update_v( 0.5 * m_dt, av[i] );
				m_Particles[i].update_x( m_dt );

			}
		}
		else{

			for( size_t i = 0; i < nParticles; i++ ){

				m_Particles[i].update_u( m_dt, sym, dudt[i], u_min[i] );

				if ( m_Settings.get<DensityCalculation>() == DensityCalculation::ContinuityEquation ){
					m_Particles[i].update_rho( m_dt, sym, drhodt[i], rho_min[i] );
				}

				m_Particles[i].update_v( m_dt, av[i], v_min[i] );
				m_Particles[i].update_x( m_dt );

			}
		}

		m_time += m_dt;
		m_CurrentTS++;

		if ( m_CurrentTS % save_step == 0 ){
			std::string Filename = "Particles";
			m_outParticles.writeOutput( Filename, m_time );
		}


	}

}

//----------------------------------------------------------------------------
// 	Subroutine to determine the right hand side of a differential
// 	equation in a single step for performing time integration
// 	In this routine and its subroutines the SPH algorithms are performed.

// 	itimestep	: Current timestep number 					[in]
//	dt			: Timestep 									[in]
//	ntotal 		: Number of particles 						[in]
//	hsml 		: Smoothing Length 							[in]
//	mass 		: Particle masses 							[in]
//	x 			: Particle position 						[in]
//	vx 			: Particle velocity 						[in]
//	u 			: Particle internal energy 					[in]
//	s 			: Particle entropy (not used here) 			[in]
//	rho 		: Density 									[in/out]
//	p 			: Pressure 									[out]
//	t 			: Temperature 								[in/out]
//	tdsdt 		: Production of viscous entropy t*ds/dt 	[out]
//	dx 			: dx = vx = dx/dt 							[out]
//	dvx 		: dvx = dvx/dt, force per unit mass 		[out]
//	du 			: du = du/dt 								[out]
//	ds 			: ds = ds/dt 								[out]
//	drho 		: drho = drho/dt 							[out]
//	itype 		: Type of particle 							[in]
//	av 			: Monaghan average velocity 				[out]
 //----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::single_step(	vNDdf<DIM>	& dvdt,
								vNdf 		& dudt,
								vNdf 		& dsdt,
								vNdf 		& tdsdt,
								vNdf		& drhodt,
								vNDdf<DIM>	& av ){

	const size_t nParticles 	= m_Particles.size();

	vNDdf<DIM> indvdt ( nParticles, vDdf<DIM>() );
	vNDdf<DIM> exdvdt ( nParticles, vDdf<DIM>() );
	vNDdf<DIM> ardvdt ( nParticles, vDdf<DIM>() );

	vNdf avdudt	(nParticles, .0 );
	vNdf ahdudt (nParticles, .0 );

	// Interaction parameters, calculating neighboring particles
	// and optimizing smoothing length
	m_pS->find();
	m_pvS->find();

	std::cout << std::endl;
	std::cout << "nInteractions:	" << m_Interactions.size() << std::endl;

	comp_density  (drhodt);

	// Internal forces:
	vNdf 	dedt( nParticles, .0 );
	// Pressure from equation of state
	for ( auto & iParticle : m_Particles ) {
		iParticle.eos();
	}

	int_force( tdsdt, dedt, indvdt );

	// Artificial
	auto artQs 	= m_Settings.getAsVector<Artificial>();
	for ( const auto iQ : artQs ){

		switch ( iQ ) {

		case Artificial::Viscosity:
			art_visc( ardvdt, avdudt);
			std::cout << "Artificial Viscosity" << std::endl;
			break;

		case Artificial::Heat:
			art_heat( ahdudt );
			std::cout << "Artificial Heat" << std::endl;
			break;

		case Artificial::NONE:
			break;

		}
	}

	// External forces:
	ext_force( exdvdt );

	// Calculating the neighboring particles and updating HSML
	if ( m_Settings.get<SmoothingLengthEvolution>() != SmoothingLengthEvolution::NONE ){
		h_upgrade( drhodt  );
	}



	// Calculating average velocity of each particle for avoiding penetration
	if ( m_Settings.get<VelocityAveraging>() != VelocityAveraging::NONE ){
		av_vel( av );
	}

	// Convert velocity, force, and energy to f and dfdt
	for( size_t i = 0; i < nParticles; i++ ){
		dvdt[i] = indvdt[i] + exdvdt[i] + ardvdt[i];
		m_Particles[i].set_acc(dvdt[i]);
		dudt[i] = dedt[i] + avdudt[i] + ahdudt[i];
	}

	if ( (m_CurrentTS+1) % m_Settings.outS().PrintStep() == 0 ){

		std::cout << "**************************************************" << std::endl;
		std::cout << "Current number of time step = " << m_CurrentTS+1 << "  	current time = " << m_time+m_dt << std::endl;
		std::cout << "**************************************************" << std::endl;

		std::cout << std::endl;
		vNsu moniParticles = m_Settings.outS().MoniParticles();
		for ( size_t i = 0; i < moniParticles.size(); i++ ){
			const size_t iMoni 	= moniParticles[i];

			std::cout << "**** Information for particle **** " << iMoni << std::endl;
			std::cout << m_Particles[iMoni];
			std::cout << "internal a:   " << indvdt[iMoni]	<< std::endl;
			std::cout << "artificial a: " << ardvdt[iMoni] 	<< std::endl;
			std::cout << "external a:   " << exdvdt[iMoni] 	<< std::endl;
			std::cout << "total a:      " << dvdt[iMoni]	<< std::endl;
			std::cout << "av:           " << av[iMoni] 		<< std::endl;
			std::cout << "dudt:         " << dudt[iMoni] 	<< std::endl;
			std::cout << "avdudt:       " << avdudt[iMoni] 	<< std::endl;
			std::cout << "ahdudt:       " << ahdudt[iMoni] 	<< std::endl;
			std::cout << "drhodt:       " << drhodt[iMoni] 	<< std::endl;
			std::cout << std::endl;

		}
	}

}

// Density approximation or change rate
template < uint8_t DIM>
void sphSim<DIM>::comp_density ( vNdf & drhodt ){

	const size_t nParticles 	= m_Particles.size();
	vNdf rho 	( nParticles, .0 );
	uint32_t i = 0;
	MYASSERT ( 	m_Settings.get<DensityCalculation>() != DensityCalculation::NONE ,
				"No density calculation specified!");

	switch ( m_Settings.get<DensityCalculation>() ) {

	case DensityCalculation::SummationDensity:
		sum_density(rho);
		i = 0;
		for ( auto & iParticle : m_Particles ) {
			iParticle.set_rho( rho[i++] );
		}
		break;

	case DensityCalculation::CSPM:
		cspm_density(rho);
		i = 0;
		for ( auto & iParticle : m_Particles ) {
			iParticle.set_rho( rho[i++] );
		}
		break;

	case DensityCalculation::ContinuityEquation:
		con_density(drhodt);
		break;

	default:
		break;
	}

}

// Summation density
//----------------------------------------------------------------------------
// 	Subroutine to calculate the density with SPH summation algorithm.

// 	ntotal 	: Number of particles 							[in]
// 	hsml	: Smoothing Length 								[in]
// 	mass 	: Particle masses 								[in]
// 	niac 	: Number of interaction pairs 					[in]
// 	pair_i 	: List of first partner of interaction pair 	[in]
// 	pair_j 	: List of second partner of interaction pair 	[in]
// 	w 		: Kernel for all interaction pairs 				[in]
// 	itype 	: type of particles 							[in]
// 	x 		: Coordinates of all particles 					[in]
// 	rho 	: Density 										[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::sum_density( vNdf & rho ){

	// Secondly calculate the rho integration over the space
	uint32_t i = 0;
	for ( auto & iParticle : m_Particles ) {
		rho[i] = iParticle.self_rho();
		i++;
	}

	// Calculate SPH sum for rho
	for ( auto & kContact : m_Interactions ) {
		kContact.calc_rho ( rho );
	}
}

template < uint8_t DIM >
void sphSim<DIM>::cspm_density( vNdf & rho ) {

	vNdf wi;

	sum_density( rho );

	// Firstly calculate the integration of the kernel over the space
	uint32_t i = 0;
	for ( auto & iParticle : m_Particles ) {
		wi[i] 	= iParticle.self_W();
		i++;
	}

	for ( auto & kContact : m_Interactions ) {
		kContact.calc_wi( wi );
	}

	// Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
	for ( size_t i = 0; m_Particles.size(); i++ ) {
		rho[i] /= wi[i];

	}

}

// Density with continuity approach
//----------------------------------------------------------------------------
// 	Subroutine to calculate 'the density with SPH continuity approach.

// 	ntotal	: Number of particles 								[in]
// 	mass	: Particle masses 									[in]
// 	niac	: Number of interaction pairs 						[in]
// 	pair_i	: List of first partner of interaction pair 		[in]
// 	pair_j	: List of second partner of interaction pair		[in]
// 	dwdx	: derivation of Kernel for all interaction pairs 	[in]
// 	vx		: Velocities of all particles 						[in]
// 	itype	: type of particles 								[in]
// 	x		: Coordinates of all particles 						[in]
// 	rho		: Density 											[in]
// 	drhodt	: Density change rate of each particle 				[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::con_density( vNdf & drhodt ) {

	size_t nParticles 	= m_Particles.size();
	drhodt = vNdf(nParticles, .0);

	for ( auto & kContact : m_Interactions ) {
		kContact.calc_con( drhodt );
	}

}


// Calculate internal forces
//----------------------------------------------------------------------------
// 	Subroutine to calculate the internal forces on the right hand side
// 	of the Navier-Stokes equations, i.e. the pressure gradient and the
// 	gradient of the viscous stress tensor, used by the time integration,
// 	Moreover the entropy production due to viscous dissipation, tds/dt,
// 	and the change of internal energy per mass, de/dt, are calculated.

// 	itimestep	: Current timestep number							[in]
//	dt			: Time step											[in]
//	ntotal 		: Number of particles								[in]
//	hsml		: Smoothing Length									[in]
//	mass		: Particle masses									[in]
//	vx			: Velocities of all particles						[in]
//	niac		: Number of interaction pairs						[in]
//	rho			: Density											[in]
//	eta			: Dynamic viscosity									[in]
//	pair_i 		: List of first partner of interaction pair			[in]
//	pair_j 		: List of second partner of interaction pair		[in]
//	dwdx		: Derivative of kernel with respect to x, y and z	[in]
//	itype 		: Type of particle (material types)					[in]
//	u			: Particle internal energy							[in]
//	x			: Particle coordinates								[in]
//	itype 		: Particle type										[in]
//	t			: Particle temperature								[in/out]
//	c			: Particle sound speed								[out]
//	p			: Particle pressure									[out]
//	dvxdt 		: Acceleration with respect to x, y and z			[out]
//	tdsdt 		: Production of viscous entropy						[out]
//	dedt		: Change of specific internal energy				[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::int_force(vNdf  		& tdsdt,
							vNdf  		& dedt,
							vNDdf<DIM> 	& dvdt ) const {


	const size_t nParticles 	= m_Particles.size();

	vNSdf<DIM> tMatrix(nParticles, symdf<DIM>() ) ;

	vNdf 	vcc(nParticles,.0);

	const bool isViscous = m_Settings.get<PDEs>() == PDEs::NavierStokesEquations;

	// Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab, vc, c
	if ( isViscous ){

		for ( auto & kInteraction : m_Interactions ) {
			kInteraction.calc_tMatrix ( tMatrix , vcc );
		}
	}

	if ( isViscous ){
		uint32_t i = 0;
		for ( auto & iParticle : m_Particles ) {
			//Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
			tdsdt[i] = iParticle.update_tdsdt( tMatrix[i] );
			i++;
		}
	}


	// Calculate SPH sum for pressure force -p,a/rho
	// and viscous force (eta Tab),b/rho
	// and the internal energy change de/dt due to -p/rho vc,c

	switch ( m_Settings.get<ParticleApproximation>() ) {

	case ParticleApproximation::Algorithm_1:
		std::cout << "::Algorithm_1" << std::endl;
		for ( auto & kInteraction : m_Interactions ) {
			kInteraction.SPH_Algorithm1( tMatrix, isViscous, dedt, dvdt );
		}
		break;

	case ParticleApproximation::Algorithm_2:
		for ( auto & kInteraction : m_Interactions ) {
			kInteraction.SPH_Algorithm2( tMatrix, isViscous, dedt, dvdt );
		}
		break;

	case ParticleApproximation::NONE:
		MYASSERT (false, "No SPH Algorithm specified!");
	}

	// Change of specific internal energy de/dt = T ds/dt - p/rho vc,c:
	for( size_t i = 0; i < nParticles; i++ ){
		dedt[i] = tdsdt[i] + 0.5 * dedt[i];
	}

}


// Calculate the artificial viscosity
//----------------------------------------------------------------------------
// 	Subroutine to calculate the artificial viscosity (Monaghan, 1992)

// 	ntotal	: Number of particles (including virtual particles)	[in]
// 	hsml	: Smoothing Length									[in]
// 	mass	: Particle masses									[in]
//	x		: Coordinates of all particles						[in]
//	vx		: Velocities of all particles						[in]
// 	niac	: Number of interaction pairs						[in]
// 	rho		: Density											[in]
// 	c		: Temperature										[in]
// 	pair_i	: List of first partner of interaction pair			[in]
// 	pair_j	: List of second partner of interaction pair		[in]
// 	w		: Kernel for all interaction pairs					[in]
// 	dwdx	: Derivative of kernel with respect to x, y and z	[in]
// 	dvxdt	: Acceleration with respect to x, y and z			[out]
// 	dedt	: Change of specific internal energy				[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::art_visc(	vNDdf<DIM> 	& dvdt,
							vNdf  		& dedt){

	const size_t nParticles = m_Particles.size();
	dvdt = vNDdf<DIM> ( nParticles, vDdf<DIM>() );
	dedt = vNdf	 ( nParticles, .0 );

	// Calculate SPH sum for artificial viscosity
	for ( auto & kInteraction : m_Interactions ) {
		kInteraction.art_viscosity( dvdt, dedt );
	}

	// Change of specific internal energy
	dedt	*= 0.5;

//	std::cout << "sphSim<DIM>::art_visc: dvdt: " <<  dvdt[0] << std::endl;
//	std::cout << "sphSim<DIM>::art_visc: dedt: " <<  dedt[0] << std::endl;

}

// Calculate the artificial heat
//----------------------------------------------------------------------------
// 	Subroutine to calculate the artificial heat, Fulk, 1994, p, a-17)

// 	ntotal	: Number of particles								[in]
// 	hsml	: Smoothing Length									[in]
// 	mass	: Particle masses									[in]
// 	x		: Coordinates of all particles						[in]
// 	vx		: Velocities of all particles						[in]
// 	rho		: Density											[in]
// 	u		: specific internal energy							[in]
// 	c		: Sound velocity									[in]
// 	niac	: Number of interaction pairs						[in]
//	pair_i	: List of first partner of interaction pair			[in]
//	pair_j	: List of second partner of interaction pair		[in]
//	w		: Kernel for all interaction pairs					[in]
//	dwdx	: Derivative of kernel with respect to x, y and z	[in]
//	dedt	: produced artificial heat, adding to energy Eq.	[out]

//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::art_heat( vNdf  & dedt ){

	const size_t nParticles = m_Particles.size();
	dedt = vNdf	 ( nParticles, .0 );
	vNdf vcc ( nParticles, .0 );

	for ( auto & kInteraction : m_Interactions ) {
		kInteraction.calc_vcc ( vcc );
	}

	for ( auto & kInteraction : m_Interactions ) {
		kInteraction.art_heat( vcc , dedt );
	}

	for( size_t i = 0; i < nParticles; i++ ){
		dedt[i] *= 2.0;
	}

}


// Calculate the external forces
//----------------------------------------------------------------------------
// 	Subroutine to calculate the external forces, e.g. gravitational forces,
//	The forces from the interactions with boundary virtual particles
//	are also calculated here as external forces.

// 	ntotal 	: Number of particles 							[in]
// 	mass	: Particle masses 								[in]
// 	x		: Coordinates of all particles					[in]
//	pair_i 	: List of first partner of interaction pair		[in]
// 	pair_j 	: List of second partner of interaction pair	[in]
// 	itype	: type of particles								[in]
// 	hsml	: Smoothing Length								[in]
//	dvxdt	: Acceleration with respect to x, y and z		[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::ext_force( vNDdf<DIM> & dvdt ){

	const size_t nParticles = m_Particles.size();

	const std::vector<ExternalForce> extForces = m_Settings.getAsVector<ExternalForce>();

	dvdt = vNDdf<DIM> ( nParticles, vDdf<DIM>() );

	for ( const auto extForce : extForces ){

		switch ( extForce ){

		case ExternalForce::Gravity: {
			vDdf<DIM> g = vDdf<DIM>();
			g[DIM] 	= -9.8;
			dvdt += vNDdf<DIM> ( nParticles, g );
			break;
		}

		case ExternalForce::BoundaryForce:
			// Boundary particle force and penalty anti-penetration force.
			for ( auto & kInteraction : m_Interactions ) {
				kInteraction.boundary_force( dvdt );
			}
			break;

		case ExternalForce::NONE:
			continue;

		}
	}

}


// Evolve the smoothing length
//----------------------------------------------------------------------------
//	Subroutine to evolve smoothing length

// 	dt		: time step											[in]
// 	ntotal 	: Number of particles 								[in]
//	mass	: Particle masses									[in]
//	vx		: Velocities of all particles						[in]
//	rho		: Density											[in]
//	niac	: Number of interaction pairs						[in]
//	pair_i	: List of first partner of interaction pair			[in]
//	pair_j 	: List of second partner of interaction pair		[in]
//	dwdx	: Derivative of kernel with respect to x, y and z	[in]
//	hsml	: Smoothing Length									[in/out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::h_upgrade( const vNdf & drhodt  ){

	size_t i = 0;
	for ( auto & iParticle : m_Particles ) {
		iParticle.update_h( drhodt[i] , m_dt );
		i++;
	}

}


//----------------------------------------------------------------------------
// 	Subroutine to calculate the average velocity to correct velocity
// 	for preventing penetration (Monaghan, 1992)

// 	ntotal	: Number of particles 							[in]
// 	mass	: Particle masses								[in]
// 	niac	: Number of interaction pairs					[in]
// 	pair_i	: List of first partner of interaction pair		[in]
// 	pair_j	: List of second partner of interaction pair	[in]
// 	w		: Kernel for all interaction pairs				[in]
// 	vx		: Velocity of each particle						[in]
// 	rho		: Density of each particle						[in]
// 	av		: Average velocity of each particle				[out]
//----------------------------------------------------------------------------
template < uint8_t DIM >
void sphSim<DIM>::av_vel( vNDdf<DIM> & av ){

	// epsilon --- a small constants chosen by experience, may lead to instability.
	// for example, for the 1 dimensional shock tube problem, the E <= 0.3

	constexpr double epsilon = 0.3;

	size_t nParticles = m_Particles.size();
	av = vNDdf<DIM>( nParticles , vDdf<DIM>() );

	for ( auto & kInteraction : m_Interactions ) {
		kInteraction.av_vel( av );
	}

	for ( size_t i = 0; i < nParticles; i++ ) {
		for ( int d = 0; d < DIM; i++){
			av[i][d] *= epsilon;
		}
	}

}

//REGISTER_3D(SPH);

template class sphSim<1>;
template class sphSim<2>;
template class sphSim<3>;

} /* namespace SPH */
































