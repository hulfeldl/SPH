//============================================================================
// Name        : time_print.cpp
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

void time_print(){

//----------------------------------------------------------------------------
//	Print out the current date and time.
//	Notes:

// 	The standard Fortran 90 routine DATE_AND_TIME is used to get
//	the current date and time strings.
//----------------------------------------------------------------------------

	using std::chrono::system_clock;

	system_clock::time_point today = system_clock::now();

	std::time_t tt;

	tt = system_clock::to_time_t ( today );
	std::cout << "Date and Time: " << ctime(&tt);

}




