//============================================================================
// Name        : main.cpp
// Author      : Lorenz
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================



#include "Settings/Settings.h"
#include "SPH/sphSim.h"


int main (int argc, char *argv[]) {

	// Create Settings
	SPH::Settings::Settings s( argc, argv );

	// Create Simulation
	switch (s.dim()) {
	case 1: {

		SPH::sphSim<1> Simulation( s );

		// Run, Forest, Run!
		Simulation.run();

		break;
	}

	case 2: {

		SPH::sphSim<2> Simulation( s );

		// Run, Forest, Run!
		Simulation.run();

		break;
	}

	case 3: {

		SPH::sphSim<3> Simulation( s );

		// Run, Forest, Run!
		Simulation.run();

		break;
	}

	default:
		MYASSERT (false, "Dimension must be 1, 2 or 3!");
	}




	return 0;
}















