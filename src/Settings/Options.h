/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * Options.h
 *
 *  Created on: Dec 9, 2014
 *      Author: kustepha
 */

#ifndef SETTINGS_OPTIONS_H_
#define SETTINGS_OPTIONS_H_

#include <string>

namespace SPH {

namespace Settings {

class Options {

private:

	std::string m_sFileName;   // name of settings file
	std::string m_iDir = "./"; // input directory
	std::string m_oDir = "./"; // output directory

	//uint64_t npc_ = 0;        // number of particles to preallocate per core

public:

	Options ();

	Options (int argc, char* argv[]);

	virtual ~Options ();

	std::string sFilename () const;

	std::string iDir () const;

	std::string oDir () const;

};

} /* namespace Settings */
} /* namespace SPH */

#endif /* SETTINGS_OPTIONS_H_ */




















