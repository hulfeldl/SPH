/*
 * Input.h
 *
 *  Created on: Aug 6, 2020
 *      Author: hulfeldl
 */




#ifndef IO_INPUT_H_
#define IO_INPUT_H_

#include <string>
#include <vector>

#include "Primitives/DataTypes/MyTypes.h"
#include "Particles/Particle.h"
#include "Output.h"
#include <vector>
#include <array>

namespace SPH {

template < uint8_t DIM, ParticleType Type >
class Input {

private:

	const std::string m_inDir;

public:

	Input();
	Input( const std::string inDir );


	virtual ~Input();

	std::vector<Particle<DIM,Type>> get_input(	const Initialization & ini,
												const ProblemType pType = ProblemType::NONE) const;

private:

	std::vector<Particle<DIM,Type>> read( 	const std::string & Filename,
											const Format FileFormat = Format::Binary ) const;

	std::vector<Particle<DIM,Type>> readBinary( const std::string & Filename ) const;

	std::vector<Particle<DIM,Type>> readTecplot( const std::string & Filename ) const;

	static std::vector<Particle<DIM,Type>> shock_tube();

	static std::vector<Particle<DIM,Type>> shear_cavity();

};

} /* namespace SPH */

#endif /* IO_INPUT_H_ */



























