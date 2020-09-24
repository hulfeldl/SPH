//============================================================================
// Name        : virt_part.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>
#include <fstream>

#include "param.inc"

void virt_part(	const int itimestep,
				const int ntotal,
				int nvirt,
				double hsml[],
				double mass[],
				double x[][],
				double vx[][],
				double rho[],
				double u[],
				double p[],
				int itype[] ){

//----------------------------------------------------------------------------
//	Subroutine to determine the information of virtual particles
//	Here only the Monaghan type virtual particles for the 2D shear
//	cavity driven problem are generated.

// 	itimestep 	: Current time step				[in]
// 	ntotal 		: Number of particles			[in]
// 	nvirt 		: Number of virtual particles	[out]
// 	hsml		: Smoothing Length				[in/out]
// 	mass		: Particle masses				[in/out]
// 	x			: Coordinates of all particles	[in/out]
// 	vx			: Velocities of all particles	[in/out]
// 	rho			: Density						[in/out]
// 	u			: internal energy				[in/out]
// 	itype		: type of particles				[in/out]
//----------------------------------------------------------------------------

	int i, j, d, im, mp;

	double x1, dx, v_inf;

	if (vp_input){

		std::ifstream xv_file, state_file, other_file;
		xv_file.open(	"../data/xv_vp.dat",	std::ios::in | std::ios::binary );
		state_file.open("../data/state_vp.dat",	std::ios::in | std::ios::binary );
		other_file.open("../data/other_vp.dat",	std::ios::in | std::ios::binary );

		// Read ntotal and dim
		int ntotal, dim;
		xv_file.read( reinterpret_cast<char*>(&ntotal),sizeof(ntotal) );
		xv_file.read( reinterpret_cast<char*>(&dim),	sizeof(dim) );

		state_file.read( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
		state_file.read( reinterpret_cast<char*>(&dim),	sizeof(dim) );

		other_file.read( reinterpret_cast<char*>(&ntotal),	sizeof(ntotal) );
		other_file.read( reinterpret_cast<char*>(&dim),	sizeof(dim) );


		for( int j = 0; j < nvirt; j++ ){

			const int i = ntotal + j;

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

		nvirt 	= 0;
		mp 		= 40;
		x1 		= 1.0e-3;
		dx 		= x1/mp;
		v_inf 	= 1.0e-3;

		// Monaghan type virtual particle on the Upper side
		for( int i = 0; i < 2*mp + 1; i++ ){
			nvirt = nvirt + 1;
			x(1, ntotal + nvirt) = i*0.5*dx;
			x(2, ntotal + nvirt) = x1;
			vx(1,ntotal + nvirt) = v_inf;
			vx(2,ntotal + nvirt) = 0.;
		}


		// Monaghan type virtual particle on the Lower side
		for( int i = 0; i < 2*mp + 1; i++ ){
			nvirt = nvirt + 1;
			x(1, ntotal + nvirt) = i*0.5*dx;
			x(2, ntotal + nvirt) = .0;
			vx(1,ntotal + nvirt) = .0;
			vx(2,ntotal + nvirt) = .0;
		}

		// Monaghan type virtual particle on the Left side
		for( int i = 0; i < 2*mp - 1; i++ ){
			nvirt = nvirt + 1;
			x(1, ntotal + nvirt) = .0;
			x(2, ntotal + nvirt) = (i+1)*0.5*dx;
			vx(1,ntotal + nvirt) = .0;
			vx(2,ntotal + nvirt) = .0;
		}

		// Monaghan type virtual particle on the Right side
		for( int i = 0; i < 2*mp - 1; i++ ){
			nvirt = nvirt + 1;
			x(1, ntotal + nvirt) = x1;
			x(2, ntotal + nvirt) = (i+1)*0.5*dx;
			vx(1,ntotal + nvirt) = .0;
			vx(2,ntotal + nvirt) = .0;
		}

		for( int i = 0; i < nvirt; i++ ){
			rho (ntotal + i) 	= 1000.0;
			mass(ntotal + i) 	= rho(ntotal + i)*dx*dx;
			p(ntotal + i) 		= .0;
			u(ntotal + i) 		= 357.1;
			itype(ntotal + i) 	= -2;
			hsml(ntotal + i) 	= dx;
		}
	}



	if ( itimestep%save_step == 0 ){

		std::ofstream xv_file, state_file, other_file;

		xv_file.open(	"../data/xv_vp.dat",	std::ios::out | std::ios::trunc | std::ios::binary	);
		state_file.open("../data/state_vp.dat",	std::ios::out | std::ios::trunc | std::ios::binary );
		other_file.open("../data/other_vp.dat",	std::ios::out | std::ios::trunc | std::ios::binary );

		if( xv_file.is_open() && state_file.is_open() && other_file.is_open() ){

			// Write ntotal and dim
			xv_file.write( 		reinterpret_cast<char*>(&nvirt),	sizeof(nvirt) );
			xv_file.write( 		reinterpret_cast<char*>(&dim),		sizeof(dim) );

			state_file.write( 	reinterpret_cast<char*>(&nvirt),	sizeof(nvirt) );
			state_file.write( 	reinterpret_cast<char*>(&dim),		sizeof(dim) );

			other_file.write( 	reinterpret_cast<char*>(&nvirt),	sizeof(nvirt) );
			other_file.write( 	reinterpret_cast<char*>(&dim),		sizeof(dim) );

			for( int i = ntotal; i < ntotal + nvirt; i++ ){

				// write  x, xv
				xv_file.write( reinterpret_cast<char*>(&i),			sizeof(i) );
				xv_file.write( reinterpret_cast<char*>(&x[i]),		sizeof(x[i]) );
				xv_file.write( reinterpret_cast<char*>(&vx[i]),		sizeof(vx[i]) );

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
	}

	if ( itimestep%print_step == 0 ){

		if (int_stat){
			std::cout << " >> Statistics: Virtual boundary particles: " << std::endl;
			std::cout << "		Number of virtual particles: " << nvirt << std::endl;
		}
	}

}




























