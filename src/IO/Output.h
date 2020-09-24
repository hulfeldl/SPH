/*
 * Output.h
 *
 *  Created on: Aug 6, 2020
 *      Author: hulfeldl
 */




#ifndef IO_OUTPUT_H_
#define IO_OUTPUT_H_

#include <string>
#include <vector>

#include "Particles/Particle.h"

namespace SPH {

template < uint8_t DIM, ParticleType Type >
class Output {

private:

	const std::vector<Particle<DIM,Type>> & m_outParticles;
	const std::string m_outDir;

	Format m_FileFormat;

public:


	Output(	const std::vector<Particle<DIM,Type>> & outParticles,
			const std::string & outDir,
			const Format FileFormat
			 );

	virtual ~Output();

	void writeOutput( 	const std::string & Filename, const double t ) const;

	void writeOutput( 	const std::string & Filename,
						const double t,
						const Format Fileformat );

private:

	void writeBinaryOutput( const std::string & Filename ) const;

	void writeTecplot( const std::string & Filename, const double simTime  ) const;

};

} /* namespace SPH */

#endif /* IO_OUTPUT_H_ */


















