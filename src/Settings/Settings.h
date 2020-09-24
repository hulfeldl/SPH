/*
 * Settings.h
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */




#ifndef SETTINGS_SETTINGS_H_
#define SETTINGS_SETTINGS_H_

#include <vector>
#include <set>

#include "Primitives/PrimitiveFactory.h"
#include "Primitives/DataTypes/MyTypes.h"
#include "Options.h"
#include "Dictionary.h"

namespace SPH {

class Initialization {

	InitialConfiguration m_IC;
	Format m_fFormat;
	std::string m_Filename;

public:

	Initialization( const Settings::Dictionary & dict );

	InitialConfiguration IC() const { return m_IC; }
	Format Fileformat() const { return m_fFormat; }
	std::string Filename() const { return m_Filename; }

};

namespace Settings {

// Control parameters for output
class OutputSettings {

	typedef 	std::vector<int>	vNsi;

private:

	PrintStats m_pStats		= PrintStats::ALL;
	uint32_t m_PrintStep	= 1;				// print_step	: Print Timestep (On Screen)
	uint32_t m_SaveStep;						// save_step	: Save Timestep (To Disk File)
	Format m_FileFormat;
	vNsu m_MoniParticles 	= {0};				// moni_particle: The particle number for information monitoring.

public:

	OutputSettings ( const Dictionary & dict );

	virtual ~OutputSettings();

	vNsu MoniParticles() const;

	uint32_t PrintStep() const;

	uint32_t SaveStep() const;

	Format Fileformat () const;

	PrintStats pStats() const;

};

class Settings {

	typedef 	std::vector<int>	vNsi;

private:

	Options m_opts;

	DictionarySource m_dictSrc;

	Dictionary m_dict;

	//----------------------------------------------------------------------------
	//	Including file for parameters and constants used
	// 	in the entire SPH software packages.
	//----------------------------------------------------------------------------

	// dim: Dimension of the problem (1, 2 or 3)
	const uint8_t DIM;

	ProblemType m_pType;

	OutputSettings m_outS;

	// maxn: Maximum number of particles
	// max_interation : Maximum number of interaction pairs
	uint32_t m_maxParticles;
	uint32_t m_maxInteractions;

	// Algorithmic Parameters
	ParticleApproximation m_paSPH;
	NearestNeighbourSearch m_nnps;
	SmoothingLengthEvolution m_sle;
	SmoothingKernel m_skf;

	// Switches for different scenarios
	//----------------------------------------------------------------------------
	DensityCalculation m_dc;
	VelocityAveraging m_va;
	PDEs m_EquationType;
//	Symmetry m_pSym 							= Symmetry::NONE;

	std::vector<ExternalForce> 	m_extForces;
	std::vector<Artificial> 	m_artQuantities;
	Initialization m_iniPhys;
	Initialization m_iniVirt;




public:

	// Std. Constructor
	Settings( int argc, char *argv[] );

	// Std. Destructor
	virtual ~Settings();

	const Options& opts() const;

	const Dictionary& dict() const;

	uint8_t dim() const;

	ProblemType pType() const;

	const OutputSettings& outS() const;

	uint32_t maxParticles() const;

	uint32_t maxInteractions() const;

	Initialization iniPhys() const { return m_iniPhys; }

	Initialization iniVirt() const { return m_iniVirt; }

	template < typename T >
	T get( ) const;

	template < typename T >
	T get( std::string name ) const;

	template < typename T >
	std::vector<T> getAsVector() const;

	template < typename T >
	std::set<T> getAsSet() const;

};

template < typename T, typename R >
const R& checkType( const R & output){
	MYASSERT( (std::is_same<T,R>::value == true) ,
			   "Types are not the same!");
	return output;
}

template < typename T >
T Settings::get( std::string name ) const {
	MYASSERT( false , "Entry not found!");
}

template < typename T >
T Settings::get() const {
	MYASSERT( false , "Entry not found!");
}

template <>
inline
ParticleApproximation Settings::get( ) const {
	return m_paSPH;
}

template <>
inline
NearestNeighbourSearch Settings::get( ) const {
	return m_nnps;
}

template <>
inline
SmoothingLengthEvolution Settings::get( ) const {
	return m_sle;
}

template <>
inline
SmoothingKernel Settings::get( ) const {
	return m_skf;
}

template <>
inline
DensityCalculation Settings::get( ) const {
	return m_dc;
}

template <>
inline
VelocityAveraging Settings::get( ) const {
	return m_va;
}

template <>
inline
PDEs Settings::get( ) const {
	return m_EquationType;
}

template < typename T >
inline
std::vector<T> Settings::getAsVector() const {
	MYASSERT( false , "Entry not found!");
}

template <>
inline
std::vector<ExternalForce> Settings::getAsVector( ) const {
	return m_extForces;
}

template <>
inline
std::vector<Artificial> Settings::getAsVector( ) const {
	return m_artQuantities;
}

template < typename T >
inline
std::set<T> Settings::getAsSet() const {
	MYASSERT( false , "Entry not found!");
}


} /* namespace Settings */
} /* namespace SPH */

#endif /* SETTINGS_SETTINGS_H_ */










