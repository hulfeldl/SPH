//============================================================================
// Name        : time_integration.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>

#include "param.inc"

void time_integration(	double x[][],
						double vx[][],
						const double mass[],
						double rho[],
						double p[],
						double u[],
						double c[],
						double s[],
						double e[],
						const int itype[],
						double hsml[],
						const int ntotal,
						const int maxtimestep,
						const double dt ){

//----------------------------------------------------------------------------

// 	x 			: coordinates of particles				[in/out]
// 	vx 			: velocities of particles				[in/out]
// 	mass 		: mass of particles						[in]
//	rho			: dnesities of particles				[in/out]
//	p			: pressure of particles					[in/out]
//	u			: internal energy of particles			[in/out]
//	c			: sound velocity of particles			[out]
//	s			: entropy of particles, not used here	[out]
//	e			: total energy of particles				[out]
//	itype		: types of particles					[in]
//					= 1 ideal gas
//					= 2 water
// 					= 3 tnt
// 	hsml		: smoothing lengths of particles		[in/out]
// 	ntotal		: total particle number					[in]
// 	maxtimestep	: maximum timesteps						[in]
//	dt			: timestep								[in]
//----------------------------------------------------------------------------


	int current_ts, nstart;

	double x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),rho_min(maxn), dx(dim,maxn), dvxfdim, maxn), du(maxn),drho(maxn), av(dim, maxn), ds(maxn),t(maxn), tdsdt(maxn);

	double time, temp_rho, temp_u;

	for( int i = 0; i < ntotal; i++ ){
		for( int d = 0; d < dim; d++ ){
			av(d,i) = .0;
		}
	}

	for( int itimestep = nstart; itimestep < nstart + maxtimestep; itimestep++ ){

		current_ts	= current_ts + 1;

		if ( itimestep%print_step == 0 ){
			std::cout << " **************************************************" << std::endl;
			std::cout <<  " Current number of time step = " << itimestep << "  	current time = " << time + dt << std::endl;
			std::cout << " **************************************************" << std::endl;
		}

		// If not first time step, then update thermal energy, density and
		// velocity half a time step

		if ( itimestep != 0 ){

			for( int i = 0; i < ntotal; i++ ){

				u_min(i) 	= u(i);
				temp_u		= .0;

				if ( dim == 1 ){
					temp_u	= -nsym*p(i)*vx(1,i)/x(1,i)/rho(i);
				}

				u(i) = u(i) + 0.5*dt*( du(i) + temp_u );

				if( u(i) < 0 ){
					u(i) = .0;
				}

				if (!summation_density){

					rho_min(i) 	= rho(i);
					temp_rho	= .0;

					if ( dim == 1 ){
						temp_rho = -nsym*rho(i)*vx(1,i)/x(1,i);
					}

					rho(i) = rho(i) + 0.5*dt*( drho(i)+ temp_rho );
				}

				for( int d = 0; d < dim; d++ ){
					v_min(d,i) 	= vx(d,i);
					vx(d,i) 	= vx(d,i) + 0.5*dt*dvx(d,i);
				}
			}
		}

		// Definition of variables out of the function vector:
		single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u, s,rho, p, t, tdsdt, dx, dvx, du, ds, drho,itype, av);

		if (itimestep == 0){

			for( int i = 0; i < ntotal; i++ ){

				temp_u = .0;

				if ( dim == 1 ){
					temp_u = -nsym*p(i)*vx(l,i)/x(l,i)/rho(i);
				}

				u(i) = u(i) + 0.5*dt*( du(i) + temp_u );

				if ( u(i) < .0 ){
					u(i) = .0;
				}

				if (!summation_density ){

					temp_rho = .0;

					if ( dim == 1 ){
						temp_rho = -nsym*rho(i)*vx(1,i)/x(1,i);
					}

					rho(i) = rho(i) + 0.5*dt*( drho(i) + temp_rho );
				}

				for( int d = 0; d < dim; d++ ){
					vx(d,i) = vx(d,i) 	+ 0.5*dt*dvx(d,i) 	+ av(d,i);
					x(d,i)	= x(d,i) 	+     dt*vx(d,i);
				}

			}
		}
		else{

			for( int i = 0; i < ntotal; i++ ){

				temp_u = .0;

				if ( dim == 1 ){
					temp_u = -nsym*p(i)*vx(1,i)/x(1,i)/rho(i);
				}

				u(i) = u_min(i) + dt*( du(i) + temp_u );

				if (u(i) < .0 ){
					u(i) = .0;
				}

				if (!summation_density){

					temp_rho=0.

					if ( dim == 1 ){
						temp_rho = -nsym*rho(i)*vx(1,i)/x(1,i);
					}

					rho(i) = rho_min(i) + dt*( drho(i) + temp_rho );
				}

				for( int d = 0; d < dim; d++ ){
					vx(d,i) = v_min(d,i) 	+ dt * dvx(d,i)	+ av(d,i);
					x(d,i) 	= x(d,i) 		+ dt * vx(d,i);
				}

			}
		}

		time = time + dt;

		if ( itimestep%save_step == 0 ){
			output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal);
		}



		if ( itimestep%print_step == 0 ){
			std::cout << std::endl;
			std::cout << "x, velocity, dvx" << std::endl;
			std::cout << x(1,moni_particle) << " " << vx(1,moni_particle) << " " << dvx(1,moni_particle) << std::endl << std::endl;
		}




	}

	nstart = current_ts;

}



