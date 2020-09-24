/*
 * Input.cpp
 *
 *  Created on: Aug 6, 2020
 *      Author: hulfeldl
 */

#include "Input.h"

namespace SPH {

template < uint8_t DIM, ParticleType Type >
Input<DIM,Type>::Input() :
	m_inDir ( "" ) {

}

template < uint8_t DIM, ParticleType Type >
Input<DIM,Type>::Input( const std::string inDir ) :
	m_inDir ( inDir ) {

}

template < uint8_t DIM, ParticleType Type >
Input<DIM,Type>::~Input() = default;


//----------------------------------------------------------------------------
//	Subroutine for loading or generating initial particle information

//	x 		: coordinates of particles 			[out]
// 	vx 		: velocities of particles			[out]
// 	mass 	: mass of particles					[out]
// 	rho 	: dnesities of particles			[out]
// 	p 		: pressure of particles				[out]
// 	u 		: internal energy of particles		[out]
// 	itype 	: types of particles				[out]
// 	hsml 	: smoothing lengths of particles	[out]
// 	ntotal 	: total particle number				[out]
//----------------------------------------------------------------------------
template < uint8_t DIM, ParticleType Type >
std::vector<Particle<DIM,Type>> Input<DIM,Type>::get_input(	const Initialization & ini,
															ProblemType pType  ) const {


	std::vector<Particle<DIM,Type>> Particles;

	// load initial particle information from external disk file
	if( ini.IC() == InitialConfiguration::Load ){

		std::string Filename;

		switch ( pType ) {

			case ProblemType::ShockTube:
			case ProblemType::ShearCavity:
				Filename 	= ( Type == ParticleType::Physical ) ? "Initial" : "Initial_vP";
				break;

			case ProblemType::NONE:
				Filename = ini.Filename();
		}

		Particles 	= read( Filename, ini.Fileformat() );

		std::cout << " **************************************************" 	<< std::endl;
		std::cout << "	Loaded initial particle configuration...         "	<< std::endl;
		std::cout << "	Total number of particles	" << Particles.size() 	<< std::endl;
		std::cout << " **************************************************" 	<< std::endl;

	}
	// Generate Initial Input
	else if ( ini.IC() == InitialConfiguration::Generate ){

		switch ( pType ) {

			case ProblemType::ShockTube:
				Particles = shock_tube();
				break;

			case ProblemType::ShearCavity:
				Particles = shear_cavity();
				break;

			case ProblemType::NONE:
				MYASSERT (false, "No Proplem Type specified! Cannot generate initial configuration!");
		}

		std::cout << " **************************************************" << std::endl;
		std::cout << "	Initial particle configuration generated		"  << std::endl;
		std::cout << "	Total number of particles	" << Particles.size()  << std::endl;
		std::cout << " **************************************************" << std::endl;

		std::string Filename 	= ( Type == ParticleType::Physical ) ? "Initial" : "Initial_vP";

		Output<DIM,Type> outParticles( Particles, m_inDir, ini.Fileformat() );
		outParticles.writeOutput( Filename, .0);

		std::cout << " **************************************************" << std::endl;
		std::cout << "	Initial particle configuration saved to " << Filename  << std::endl;
		std::cout << " **************************************************" << std::endl;

	}

	return Particles;

}

template < uint8_t DIM, ParticleType Type >
std::vector<Particle<DIM,Type>> Input<DIM,Type>::read( 	const std::string & Filename,
														const Format FileFormat ) const {

	switch ( FileFormat ){

	case Format::Binary:
		return readBinary( Filename );
		break;

	case Format::Tecplot:
		 return readTecplot( Filename );
		 break;

	case Format::NONE:
		MYASSERT (false, "No Fileformat specified! Cannot load file!");
	}
}


template < uint8_t DIM, ParticleType Type >
std::vector<Particle<DIM,Type>> Input<DIM,Type>::readBinary( const std::string & Filename ) const {

	std::string FileExt	= ".bin";

	std::ifstream xv_file, state_file, other_file;

	xv_file.open(	m_inDir + Filename + "_xv"		+ FileExt,	std::ios::in | std::ios::binary );
	state_file.open(m_inDir + Filename + "_state"	+ FileExt,	std::ios::in | std::ios::binary );
	other_file.open(m_inDir + Filename + "_other"	+ FileExt,	std::ios::in | std::ios::binary );

	// Read ntotal and dim
	int ntotal, dim;
	xv_file.read ( reinterpret_cast<char*>(&ntotal),sizeof(ntotal) );
	xv_file.read ( reinterpret_cast<char*>(&dim),	sizeof(dim) );
	MYASSERT ( 	DIM == static_cast<uint8_t>(dim),
				"Dimension of input file does not match dimension of Simulation!");

	state_file.read ( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
	state_file.read ( reinterpret_cast<char*>(&dim),	sizeof(dim) );
	MYASSERT ( 	DIM == static_cast<uint8_t>(dim),
				"Dimension of input file does not match dimension of Simulation!");

	other_file.read ( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
	other_file.read ( reinterpret_cast<char*>(&dim),	sizeof(dim) );
	MYASSERT ( 	DIM == static_cast<uint8_t>(dim),
				"Dimension of input file does not match dimension of Simulation!");

	std::vector<Particle<DIM,Type>> Particles;
	for( int i = 0; i < ntotal; i++){

		int im; // dummy var for i

		// write  x, xv
		std::array<double,DIM> x;
		std::array<double,DIM> v;
		xv_file.read( reinterpret_cast<char*>( &im ),	sizeof( im ) );
		xv_file.read( reinterpret_cast<char*>( &x ),	sizeof( x ) );
		xv_file.read( reinterpret_cast<char*>( &v ),	sizeof( v ) );

		// write mass, rho, p, u
		double rho;
		double m;
		double p;
		double u;
		state_file.read( reinterpret_cast<char*>( &im ),	sizeof( im ) );
		state_file.read( reinterpret_cast<char*>( &m ),		sizeof( m ) );
		state_file.read( reinterpret_cast<char*>( &rho ),	sizeof( rho ) );
		state_file.read( reinterpret_cast<char*>( &p ),		sizeof( p ) );
		state_file.read( reinterpret_cast<char*>( &u ),		sizeof( u ) );
		const double c					= .0;
		const double s					= .0;
		const double e					= .0;

		// write itype, hsml
		int it;
		double hsml;
		state_file.read( reinterpret_cast<char*>( &im ),	sizeof( im ) );
		other_file.read( reinterpret_cast<char*>( &it ),	sizeof( it ) );
		other_file.read( reinterpret_cast<char*>( &hsml ),	sizeof( hsml ) );

		const FluidType ftype 			= static_cast<FluidType> ( it );
		const SmoothingKernel kType 	= SmoothingKernel::CubicSpline;

		const Particle<DIM,Type> iParticle( vDdf<DIM>(x), vDdf<DIM>(v), m, rho, p, u, c, s, e, kType, hsml, ftype);
		Particles.push_back( iParticle );
	}

	xv_file.close();
	state_file.close();
	other_file.close();

	return Particles;

}

template < uint8_t DIM, ParticleType Type >
std::vector<Particle<DIM,Type>> Input<DIM,Type>::readTecplot( const std::string & Filename ) const {

	std::string FileExt	= ".plt";

	std::ifstream pltFile;

	std::vector<Particle<DIM,Type>> Particles;

	return Particles;
}

template < uint8_t DIM, ParticleType Type >
std::vector<Particle<DIM,Type>> Input<DIM,Type>::shock_tube() {
	std::string pType = ( Type == ParticleType::Physical ) ? "physical" : "virtual";
	MYASSERT ( 	false,
				"Generation of " + pType + " particles initial configuration for the " +
				std::to_string(DIM) + "D shock tube problem failed! Not implemented yet!" );
}


// 1D shock tube
//----------------------------------------------------------------------------
// 	This subroutine is used to generate initial data for the
// 	I d noh shock tube problem

// 	x		: coordinates of particles			[out]
// 	vx		: velocities of particles			[out]
//	mass	: mass of particles					[out]
//	rho		: dnesities of particles			[out]
//	p		: pressure of particles				[out]
//	u		: internal energy of particles		[out]
//	itype	: types of particles				[out]
//				= 1, ideal gas
//	hsml	: smoothing lengths of particles	[out]
//	ntotal	: total particle number				[out]
//----------------------------------------------------------------------------
template <>
std::vector<Particle<1,ParticleType::Physical>> Input<1,ParticleType::Physical>::shock_tube() {

	const int ntotal		= 400;
	const double space_x	= 0.6 / 80.0;

	std::vector<Particle<1,ParticleType::Physical>> Particles;

	for( int i = 0; i < ntotal; i++ ){

		const double m					= 0.75/400.0;
		const double c					= .0;
		const double s					= .0;
		const double e					= .0;

		const SmoothingKernel kType 	= SmoothingKernel::CubicSpline;
		const double hsml 				= 0.015;

		const FluidType itype 			= FluidType::IdealGas;

		const vDdf<1> v(.0);

		if ( i < 320 ){

			const vDdf<1> x( -0.6 + space_x / 4.0 * static_cast<double>(i) );

			const double u					= 2.5;
			const double rho				= 1.0;
			const double p					= 1.0;

			const Particle<1,ParticleType::Physical> iParticle( x, v, m, rho, p, u, c, s, e, kType, hsml, itype );
			Particles.push_back( iParticle );
		}
		else{

			const vDdf<1> x( 0.0 + space_x*( static_cast<double>(i) - 319 ) );

			const double u					= 1.795;
			const double rho				= 0.25;
			const double p					= 0.1795;

			const Particle<1,ParticleType::Physical> iParticle( x, v, m, rho, p, u, c, s, e, kType, hsml, itype );
			Particles.push_back( iParticle );
//			std::cout << std::endl;

		}
	}

	return Particles;

}

template < uint8_t DIM, ParticleType Type >
std::vector<Particle<DIM,Type>> Input<DIM,Type>::shear_cavity() {
	std::string pType = ( Type == ParticleType::Physical ) ? "physical" : "virtual";
	MYASSERT ( 	false,
				"Generation of " + pType + " particles initial configuration for the " +
				std::to_string(DIM) + "D shear_cavity problem failed! Not implemented yet!" );
}

// 2D shear cavity physical particles
//----------------------------------------------------------------------------
//	This subroutine is used to generate initial data for the
//	2 d shear driven cavity probem with Re = 1

//	x		: coordinates of particles			[out]
//	vx		: velocities of particles			[out]
//	mass	: mass of particles					[out]
//	rho		: dnesities of particles			[out]
//	p		: pressure of particles				[out]
//	u		: internal energy of particles		[out]
//	itype	: types of particles				[out]
//				= 2, water
//	h 		: smoothing lengths of particles	[out]
//	ntotal	: total particle number				[out]
//----------------------------------------------------------------------------
template <>
std::vector<Particle<2,ParticleType::Physical>> Input<2,ParticleType::Physical>::shear_cavity() {


	// Giving mass and smoothing length as well as other data.
	const int m 		= 41;
	const int n 		= 41;
	const int mp 		= m - 1;
	const int np 		= n - 1;

	const double x1		= 1.0e-3;
	const double y1		= 1.0e-3;
	const double dx 	= x1/static_cast<double>(mp);
	const double dy 	= y1/static_cast<double>(np);


	std::vector<Particle<2,ParticleType::Physical>> Particles;

	for( int i = 0; i < mp; i++){
		for(int j = 0; j < np; j++){

			const vDdf<2> x( i*dx + 0.5*dx, j*dy + 0.5*dy );
			const vDdf<2> v;

			const double rho 	= 1000.0;
			const double mass	= dx*dy*rho;
			const double p		= .0;
			const double u 		= 357.1;
			const double c		= .0;
			const double s		= .0;
			const double e		= .0;

			const SmoothingKernel kType 	= SmoothingKernel::CubicSpline;
			const double hsml 				= dx;

			const FluidType itype 			= FluidType::Water;

			const Particle<2,ParticleType::Physical> iParticle( x, v, mass, rho, p, u, c, s, e, kType, hsml, itype );
			Particles.push_back( iParticle );
		}
	}

	return Particles;

}


// 2D shear cavity virtual particles
template <>
std::vector<Particle<2,ParticleType::Virtual>> Input<2,ParticleType::Virtual>::shear_cavity() {

	uint32_t nvirt 		= 0;
	const uint32_t mp 	= 40;
	const double x1 	= 1.0e-3;
	const double dx 	= x1 / static_cast<double>(mp);
	const double v_inf 	= 1.0e-3;

	vNDdf<2> x;

	// Monaghan type virtual particle on the Upper side
	for( uint32_t i = 0; i < 2*mp + 1; i++ ){
		x[nvirt] = vDdf<2>(i*0.5*dx, x1);
		nvirt++;
	}


	// Monaghan type virtual particle on the Lower side
	for( uint32_t i = 0; i < 2*mp + 1; i++ ){
		x[nvirt] = vDdf<2>( i*0.5*dx, .0);
		nvirt++;
	}

	// Monaghan type virtual particle on the Left side
	for( uint32_t i = 0; i < 2*mp - 1; i++ ){
		x[nvirt] = vDdf<2>( .0, (i+1)*0.5*dx );
		nvirt++;
	}

	// Monaghan type virtual particle on the Right side
	for( uint32_t i = 0; i < 2*mp - 1; i++ ){
		x[nvirt] = vDdf<2>( x1, (i+1)*0.5*dx );
		nvirt++;
	}

	std::vector<Particle<2,ParticleType::Virtual>> Particles;

	for( uint32_t i = 0; i < nvirt; i++ ){

		vDdf<2> v = ( x[i][0] == x1 ) ? vDdf<2>(v_inf,.0) : vDdf<2>();
		const double rho 	= 1000.0;
		const double mass	= dx*dx*rho;
		const double p		= .0;
		const double u 		= 357.1;
		const double c		= .0;
		const double s		= .0;
		const double e		= .0;

		const SmoothingKernel kType 	= SmoothingKernel::CubicSpline;
		const double hsml 				= dx;

		const FluidType itype 			= FluidType::Water;

		const Particle<2,ParticleType::Virtual> iParticle( x[i], v, mass, rho, p, u, c, s, e, kType, hsml, itype );
		Particles.push_back( iParticle );
	}

	return Particles;
}

template class Input <1,ParticleType::Physical>;
template class Input <2,ParticleType::Physical>;
template class Input <3,ParticleType::Physical>;

template class Input <1,ParticleType::Virtual>;
template class Input <2,ParticleType::Virtual>;
template class Input <3,ParticleType::Virtual>;

} /* namespace SPH */
















