/*
 * Settings.cpp
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */

#include "Settings.h"
#include <type_traits>

namespace SPH {

Initialization::Initialization( const Settings::Dictionary & dict ) :
	m_IC ( 		dict.get<InitialConfiguration>("type") ),
	m_fFormat(	dict.get<Format>("file format") ),
	m_Filename(	dict.get<std::string>("file name") ) {
}

namespace Settings {

using namespace std;

OutputSettings::OutputSettings ( const Dictionary & dict ) :
	m_pStats(		dict.get<PrintStats>("print") ),
	m_PrintStep(	dict.get<uint32_t>("step") ),
	m_SaveStep( 	dict.get<uint32_t>("save interval")),
	m_FileFormat(	dict.get<Format>("file format") ),
	m_MoniParticles(dict.getAsVector<uint32_t>("monitor particles")) {

}

OutputSettings::~OutputSettings() = default;

vNsu OutputSettings::MoniParticles() const {
	return m_MoniParticles;
}

uint32_t OutputSettings::PrintStep() const {
	return m_PrintStep;
}

uint32_t OutputSettings::SaveStep() const {
	return m_SaveStep;
}

Format OutputSettings::Fileformat () const {
	return m_FileFormat;
}

PrintStats OutputSettings::pStats() const {
	return m_pStats;
}

// Settings class
// ********************************************************************************
Settings::Settings( int argc, char *argv[] ) :
	m_opts( argc , argv ),
	m_dictSrc( m_opts.sFilename() ),
	m_dict( m_dictSrc ),
	DIM( m_dict.get<Dictionary>("simulation").get<uint8_t>("dimension") ),
	m_pType ( m_dict.get<Dictionary>("simulation").get<ProblemType>("problem type") ),
	m_outS( 	m_dict.get<Dictionary>("output").get<Dictionary>("SPH") ),
	m_maxParticles ( m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<uint32_t>("Maximal Particles") ),
	m_maxInteractions ( m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<uint32_t>("Max. Int. per Particle") ),
	m_paSPH ( 	m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<ParticleApproximation>("Particle Approximation") ),
	m_nnps  ( 	m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<NearestNeighbourSearch>("Nearest Neighbour Search") ),
	m_sle ( 	m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<SmoothingLengthEvolution>("Smoothing Length Evolution") ),
	m_skf ( 	m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<SmoothingKernel>("Smoothing Kernel") ),
	m_dc ( 		m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<DensityCalculation>("Density Calculation") ),
	m_va (		m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<VelocityAveraging>("Velocity Averaging") ),
	m_EquationType ( m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").get<PDEs>("PDEs") ),
	m_extForces ( 	m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").getAsVector<ExternalForce>("External Forces")	),
	m_artQuantities ( m_dict.get<Dictionary>("algorithm").get<Dictionary>("SPH").getAsVector<Artificial>("Artificial") ),
	m_iniPhys(	m_dict.get<Dictionary>("initialization").get<Dictionary>("physical") ),
	m_iniVirt(	m_dict.get<Dictionary>("initialization").get<Dictionary>("virtual") )
//	m_pSym (),
		{

}

Settings::~Settings() = default;

const Options& Settings::opts() const {

	return m_opts;
}

const Dictionary& Settings::dict() const {
	return m_dict;
}

uint8_t Settings::dim() const {
	return DIM;
}

ProblemType Settings::pType() const {
	return m_pType;
}

const OutputSettings& Settings::outS() const {
	return m_outS;
}

uint32_t Settings::maxParticles() const {
	return m_maxParticles;
}

uint32_t Settings::maxInteractions() const{
	return m_maxInteractions;
}



} /* namespace Settings */
} /* namespace SPH */






















