//============================================================================
// Name        : single_step.cpp
// Author      : Lorenz
// Created on  : Jul 28, 2020
// Version     :
// Copyright   : Your copyright notice
// Description : main file for SPH simulation
//============================================================================


#include <iostream>

#include "param.inc"

void single_step(	const int itimestep,
					const double dt,
					const int ntotal,
					const double hsml[],
					const double mass[],
					const double x[][],
					const double vx[][],
					const double u[],
					const double s[],
					double rho[],
					double p[],
					double t[],
					double tdsdt[],
					double dx[][],
					double dvx[][],
					double du[],
					double ds[],
					double drho[],
					const int itype[],
					double av[][] ){

//----------------------------------------------------------------------------
// 	Subroutine to determine the right hand side of a differential
// 	equation in a single step for performing time integration
// 	In this routine and its subroutines the SPH algorithms are performed.

// 	itimestep	: Current timestep number 					[in]
//	dt			: Timestep 									[in]
//	ntotal 		: Number of particles 						[in]
//	hsml 		: Smoothing Length 							[in]
//	mass 		: Particle masses 							[in]
//	x 			: Particle position 						[in]
//	vx 			: Particle velocity 						[in]
//	u 			: Particle internal energy 					[in]
//	s 			: Particle entropy (not used here) 			[in]
//	rho 		: Density 									[in/out]
//	p 			: Pressure 									[out]
//	t 			: Temperature 								[in/out]
//	tdsdt 		: Production of viscous entropy t*ds/dt 	[out]
//	dx 			: dx = vx = dx/dt 							[out]
//	dvx 		: dvx = dvx/dt, force per unit mass 		[out]
//	du 			: du = du/dt 								[out]
//	ds 			: ds = ds/dt 								[out]
//	drho 		: drho = drh,o/dt 							[out]
//	itype 		: Type of particle 							[in]
//	av 			: Monaghan average velocity 				[out]
 //----------------------------------------------------------------------------

	int nvirt, niac, pair_i(max_interaction) ,pair_j(max_interaction), ns(maxn);

	double w(max_interaction), dwdx(dim,max_interaction),indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn);

	for( int i = 0; i < ntotal; i++ ){
		avdudt(i) = .0;
		ahdudt(i) = .0;

		for( int d = 0; d < dim; d++ ){
			indvxdt(d,i) = .0;
			ardvxdt(d,i) = .0;
			exdvxdt(d,i) = .0;
		}
	}


	// Positions of virtual (boundary) particles:
	nvirt = 0;
	if ( virtual_part ){
		virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,rho,u,p,itype);
	}

	// Interaction parameters, calculating neighboring particles
	// and optimizing smoothing length
	if ( nnps == 1 ){
		direct_find(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,pair_j,w,dwdx,ns);
	}
	else if ( nnps == 2 ){
		link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,pair_j,w,dwdx,ns);
	}
	else if ( nnps == 3 ){
		tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,pair_j,w,dwdx,ns);
	}


	// Density approximation or change rate
	if (summation_density){
		sum_density(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,itype,rho);
	}
	else{
		con_density(ntotal+nvirt,mass,niac,pair_i,pair_j,dwdx,vx, itype,x,rho, drho);
	}

	// Dynamic viscosity:
	if (visc){
		viscosity(ntotal+nvirt,itype,x,rho,eta);
	}

	// Internal forces:
	int_force(itimestep,dt,ntotal+nvirt,hsml,mass,vx,niac,rho,eta, pair_i,pair_j,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du);

	// Artificial viscosity:
	if (visc_artificial){
		art_visc(ntotal+nvirt,hsml,mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt);
	}

	// External forces:
	if (ex_force){
		ext_force (ntotal+nvirt,mass,x,niac,pair_i,pair_j,itype, hsml, exdvxdt);
	}

	// Calculating the neighboring particles and updating HSML
	if (sle != O){
		h_upgrade(dt,ntotal, mass, vx, rho, niac,pair_i, pair_j, dwdx, hsml);
	}

	if (heat_artificial){
		art_heat(ntotal+nvirt,hsml,mass,x,vx,niac,rho,u, c,pair_i,pair_j,w,dwdx,ahdudt);
	}

	// Calculating average velocity of each particle for avoiding penetration
	if (average_velocity){
		av_vel(ntotal,mass,niac,pair_i,pair_j, w, vx, rho, av);
	}

	// Convert velocity, force, and energy to f and dfdt
	for( int i = 0; i < ntotal; i++ ){
		for( int d = 0; d < dim; d++ ){
			dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i);
		}
		du(i) = du(i) + avdudt(i) + ahdudt (i);
	}


	if ( itimestep%print_step == 0 ){
		std::cout <<  "**** Information for particle **** " << moni_particle << std::endl;
		std::cout <<  "internal a, artificial a, external a, total a "  << std::endl;
		std::cout <<  indvxdt(0,moni_particle) << " " << ardvxdt(0,moni_particle) << " " << exdvxdt(0,moni_particle) << " " << dvx(0,moni_particle)  << std::endl;
	}

}











