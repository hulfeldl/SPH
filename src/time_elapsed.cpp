//============================================================================
// Name        : time_elapsed.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

void time_elapsed(	std::time_t  & s){

//----------------------------------------------------------------------------
//	The standard Fortran 90 routine RTC is used to calculate the elapsed CPU
//	use dfport
//----------------------------------------------------------------------------

	using std::chrono::system_clock;

	system_clock::time_point now = system_clock::now();

	s = system_clock::to_time_t ( now );

}




