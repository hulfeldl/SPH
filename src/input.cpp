//============================================================================
// Name        : input.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>
#include <fstream>

#include "param.inc"

void input(	double x[][],
			double vx[][],
			double mass[],
			double rho[],
			double p[],
			double u[],
			int itype[],
			double hsml[],
			int & ntotal ){

//----------------------------------------------------------------------------
//	Subroutine for loading or generating initial particle information

//	x 		: coordinates of particles 			[out]
// 	vx 		: velocities of particles			[out]
// 	mass 	: mass of particles					[out]
// 	rho 	: dnesities of particles			[out]
// 	p 		: pressure of particles				[out]
// 	u 		: internal energy of particles		[out]
// 	itype 	: types of particles				[out]
// 	hsml 	: smoothing lengths of particles	[out]
// 	ntotal 	: total particle number				[out]
//----------------------------------------------------------------------------


	// load initial particle information from external disk file
	if(config_input){

		std::ifstream xv_file, state_file, other_file;
		xv_file.open(	"../data/f_xv.dat",		ios::in | ios::binary );
		state_file.open("../data/f_state.dat",	ios::in | ios::binary );
		other_file.open("../data/f_other.dat",	ios::in | ios::binary );

		std::cout << " **************************************************" << endl;
		std::cout << "	Loading initial particle configuration...         ";
		read (1,*) ntotal
		std::cout << "	Total number of particles	" << ntotal;
		std::cout << " **************************************************" << endl;

		// Read ntotal and dim
		xv_file.write( reinterpret_cast<char*>(&ntotal),sizeof(ntotal) );
		xv_file.write( reinterpret_cast<char*>(&dim),	sizeof(dim) );

		state_file.write( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
		state_file.write( reinterpret_cast<char*>(&dim),	sizeof(dim) );

		other_file.write( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
		other_file.write( reinterpret_cast<char*>(&dim),	sizeof(dim) );

		for( int i = 0; i < ntotal; i++){

			int im; // dummy var for i

			// write  x, xv
			xv_file.read( reinterpret_cast<char*>(&im),		sizeof(im) );
			xv_file.read( reinterpret_cast<char*>(&x[i]),	sizeof(x[i]) );
			xv_file.read( reinterpret_cast<char*>(&vx[i]),	sizeof(vx[i]) );

			// write mass, rho, p, u
			state_file.read( reinterpret_cast<char*>(&im),		sizeof(im) );
			state_file.read( reinterpret_cast<char*>(&mass[i]),	sizeof(mass[i]) );
			state_file.read( reinterpret_cast<char*>(&rho[i]),	sizeof(rho[i]) );
			state_file.read( reinterpret_cast<char*>(&p[i]),	sizeof(p[i]) );
			state_file.read( reinterpret_cast<char*>(&u[i]),	sizeof(u[i]) );

			// write itype, hsml
			state_file.read( reinterpret_cast<char*>(&im),		sizeof(im) );
			other_file.read( reinterpret_cast<char*>(&itype[i]),sizeof(itype[i]) );
			other_file.read( reinterpret_cast<char*>(&hsml[i]),	sizeof(hsml[i]) );
		}

		xv_file.close();
		state_file.close();
		other_file.close();

	}
	else{

		std::ofstream xv_file, state_file, other_file;

		xv_file.open(	"../data/ini_xv.dat",	ios::out | ios::trunc | ios::binary	);
		state_file.open("../data/ini_state.dat",ios::out | ios::trunc | ios::binary );
		other_file.open("../data/ini_other.dat",ios::out | ios::trunc | ios::binary );

		if (shocktube){
			shock_tube(x, vx, mass, rho, p, u,itype, hsml, ntotal);
		}

		if (shearcavity){
			shear_cavity(x, vx, mass, rho, p, u,itype, hsml, ntotal);
		}

		if( xv_file.is_open() && state_file.is_open() && other_file.is_open() ){

			// Write ntotal and dim
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
				state_file.write( reinterpret_cast<char*>(&i),	sizeof(i) );
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

		std::cout << " **************************************************" << std::endl;
		std::cout << "	Initial particle configuration generated	"  << std::endl;
		std::cout << "	Total number of particles	" << ntotal  << std::endl;
		std::cout << " **************************************************"  << std::endl;

	}

}


void shock_tube(double x[][],
				double vx[][],
				double mass[],
				double rho[],
				double p[],
				double u[],
				int itype[],
				double hsml[],
				int ntotal ){

//----------------------------------------------------------------------------
// 	This subroutine is used to generate initial data for the
// 	I d noh shock tube problem

// 	x		: coordinates of particles			[out]
// 	vx		: velocities of particles			[out]
//	mass	: mass of particles					[out]
//	rho		: dnesities of particles			[out]
//	p		: pressure of particles				[out]
//	u		: internal energy of particles		[out]
//	itype	: types of particles				[out]
//				= 1, ideal gas
//	hsml	: smoothing lengths of particles	[out]
//	ntotal	: total particle number				[out]
//----------------------------------------------------------------------------

	double space_x;

	ntotal	= 400;
	space_x	= 0.6/80.0;

	for( int i = 0; i < ntotal; i++ ){

		mass(i)		= 0.75/400.0;
		hsml(i)		= 0.015;
		itype(i)	= 1;

		for( int d = 0; d < dim; d++ ){
			x(d,i) 	= 0.0;
			vx(d,i) = 0.0;
		}
	}

	for( int i = 0; i < 320; i++ ){
		x(0,i)	= -0.6 + space_x/4.0*i;
	}

	for( int i = 320; i < ntotal; i++ ){
		x(0,i)	= 0.0 + space_x*(i-319);
	}

	for( int i = 0; i < ntotal; i++ ){
		if (x(0,i) <= 1.e-8){
			u(i)	= 2.5;
			rho(i)	= 1.0;
			p(i)	= 1.0;
		}
		else{
			u(i)	= 1.795;
			rho(i)	= 0.25;
			p(i)	= 0.1795;
		}
	}

}

void shear_cavity(	double x[][],
					double vx[][],
					double mass[],
					double rho[],
					double p[],
					double u[],
					int itype[],
					double hsml[],
					int ntotal ){

//----------------------------------------------------------------------------
//	This subroutine is used to generate initial data for the
//	2 d shear driven cavity probem with Re = 1

//	x		: coordinates of particles			[out]
//	vx		: velocities of particles			[out]
//	mass	: mass of particles					[out]
//	rho		: dnesities of particles			[out]
//	p		: pressure of particles				[out]
//	u		: internal energy of particles		[out]
//	itype	: types of particles				[out]
//				= 2, water
//	h 		: smoothing lengths of particles	[out]
//	ntotal	: total particle number				[out]
//----------------------------------------------------------------------------

	// Giving mass and smoothing length as well as other data.
	int m 	= 41;
	int n 	= 41;
	int mp 	= m-1;
	int np 	= n-1;
	ntotal	= mp*np;

	double xl	= 1.0e-3;
	double yl	= 1.0e-3;
	double dx 	= xl/mp;
	double dy 	= yl/np;

	for( int i = 0; i < mp; i++){
		for(int j = 0; j < np; j++){
			int k 	= j + i*np;
			x(0, k) = i*dx + 0.5*dx;
			x(1, k) = j*dy + 0.5*dy;
		}
	}

	for( int i = 0; i < ntotal; i++){
		vx(0, i) 	= .0;
		vx(1, i)	= .0;
		rho (i) 	= 1000.0;
		mass(i) 	= dx*dy*rho(i);
		p(i)		= .0;
		u(i) 		= 357.1;
		itype(i) 	= 2;
		hsml(i) 	= dx;
	}

}








