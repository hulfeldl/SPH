//============================================================================
// Name        : output.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>
#include <fstream>

#include "param.inc"

void output(const double x[][],
			const double vx[][],
			const double mass[],
			const double rho[],
			const double p[],
			const double u[],
			const double c[],
			const int itype[],
			const double hsml[],
			const int ntotal ){

//----------------------------------------------------------------------------
//	Subroutine for saving particle information to external disk file

// 	x 		: coordinates of particles			[in]
// 	vx		: velocities of particles			[in]
// 	mass	: mass of particles					[in]
// 	rho 	: dnesities of particles			[in]
//	p 		:pressure of particles				[in]
//	u 		: internal energy of particles		[in]
//	c 		: sound velocity of particles		[in]
//	itype 	: types of particles				[in]
//	hsml 	: smoothing lengths of particles	[in]
//	ntotal 	: total particle number				[in]
//----------------------------------------------------------------------------

	std::ofstream xv_file, state_file, other_file;

	xv_file.open(	"../data/f_xv.dat",		ios::out | ios::trunc | ios::binary	 );
	state_file.open("../data/f_state.dat",	ios::out | ios::trunc | ios::binary	 );
	other_file.open("../data/f_other.dat",	ios::out | ios::trunc | ios::binary	 );

	xv_file.write( reinterpret_cast<char*>(&ntotal),sizeof(ntotal) );
	xv_file.write( reinterpret_cast<char*>(&dim),	sizeof(dim) );

	state_file.write( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
	state_file.write( reinterpret_cast<char*>(&dim),	sizeof(dim) );

	other_file.write( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
	other_file.write( reinterpret_cast<char*>(&dim),	sizeof(dim) );

	for( int i = 0; i < ntotal; i++){

		// write  x, xv
		xv_file.write( reinterpret_cast<char*>(&i),		sizeof(i) );
		xv_file.write( reinterpret_cast<char*>(&x[i]),	sizeof(x[i]) );
		xv_file.write( reinterpret_cast<char*>(&vx[i]),	sizeof(vx[i]) );

		// write mass, rho, p, u
		state_file.write( reinterpret_cast<char*>(&i),		sizeof(i) );
		state_file.write( reinterpret_cast<char*>(&mass[i]),sizeof(mass[i]) );
		state_file.write( reinterpret_cast<char*>(&rho[i]),	sizeof(rho[i]) );
		state_file.write( reinterpret_cast<char*>(&p[i]),	sizeof(p[i]) );
		state_file.write( reinterpret_cast<char*>(&u[i]),	sizeof(u[i]) );

		// write itype, hsml
		state_file.write( reinterpret_cast<char*>(&i),			sizeof(i) );
		other_file.write( reinterpret_cast<char*>(&itype[i]),	sizeof(itype[i]) );
		other_file.write( reinterpret_cast<char*>(&hsml[i]),	sizeof(hsml[i]) );
	}

	xv_file.close();
	state_file.close();
	other_file.close();

}




