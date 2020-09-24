/*
 * Timer.cpp
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */

#include "Timer.h"

namespace SPH {

	Timer::Timer() {
		// TODO Auto-generated constructor stub

	}

	Timer::~Timer() {
		// TODO Auto-generated destructor stub
	}

	//----------------------------------------------------------------------------
	//	The standard Fortran 90 routine RTC is used to calculate the elapsed CPU
	//	use dfport
	//----------------------------------------------------------------------------
	const std::time_t Timer::time_elapsed(	){


		using std::chrono::system_clock;

		system_clock::time_point now = system_clock::now();

		std::time_t s;

		return s = system_clock::to_time_t ( now );

	}

	//----------------------------------------------------------------------------
	//	Print out the current date and time.
	//	Notes:

	// 	The standard Fortran 90 routine DATE_AND_TIME is used to get
	//	the current date and time strings.
	//----------------------------------------------------------------------------
	void Timer::time_print(){


		using std::chrono::system_clock;

		system_clock::time_point today = system_clock::now();

		std::time_t tt;

		tt = system_clock::to_time_t ( today );
		std::cout << "Date and Time: " << ctime(&tt) << std::endl;

	}

} /* namespace SPH */








