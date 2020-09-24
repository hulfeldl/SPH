/*
 * Timer.h
 *
 *  Created on: Aug 4, 2020
 *      Author: hulfeldl
 */




#ifndef UTILITIES_TIMER_H_
#define UTILITIES_TIMER_H_

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

namespace SPH {

class Timer {

public:

	Timer();

	virtual ~Timer();

	// time Elapsed
	const std::time_t time_elapsed(	);

	//	Print out the current date and time.
	void time_print();

};

} /* namespace SPH */

#endif /* UTILITIES_TIMER_H_ */










