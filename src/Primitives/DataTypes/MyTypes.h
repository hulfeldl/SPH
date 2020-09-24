/*
 * MyTypes.h
 *
 *  Created on: Aug 13, 2020
 *      Author: hulfeldl
 */




#ifndef PRIMITIVES_DATATYPES_MYTYPES_H_
#define PRIMITIVES_DATATYPES_MYTYPES_H_

#include <array>
#include <vector>

#include "Primitives/MyError.h"
#include "MyMatrix.h"
#include "MyVector.h"
#include "Symmetric.h"

namespace SPH {

typedef unsigned char     		uint8_t;

// Typedefs
typedef std::array<double,2> 	v2df;
typedef std::array<double,4> 	v4df;
typedef std::array<double,6> 	v6df;
typedef std::vector<double> 	vNdf;

typedef std::vector<int>		vNsi;
typedef std::vector<uint32_t> 	vNsu;

template < uint8_t DIM >
using vDdf = MyVector<double,DIM>;

template < uint8_t DIM >
using symdf = Symmetric<double,DIM>;

typedef std::vector<v4df> 		vN4df;

template < uint8_t DIM >
using vNDdf = std::vector<vDdf<DIM>>;

template < uint8_t DIM >
using vNSdf = std::vector<symdf<DIM>>;

#define REGISTER_3D(template_name)	template class template_name##"<1>" \
									template class template_name##"<2>" \
									template class template_name##"<3>"

//	Enums
// -------------------------------------------------------------------------------

// IO
enum class Format {
	NONE,				//	0	: File Format not specified
	Binary,				//	1	: Binary File
	Tecplot				//	2	: Tecplot File Format
};

// Particles
enum class FluidType {
	NONE, 							// 	0	: No Fluid Type Specified
	IdealGas, 						// 	1	: Ideal Gas
	Water 							// 	2 	: Water
};

// SPH parameters
// -------------------------------------------------------------------------------

// Print statistics about particle interactions
enum class PrintStats {
	NONE,				//	0 : Print nothing
	sphParticles,		// 	1 : Print only information about physical Particles
	virtParticles, 		// 	2 : Print only information about virtual Particles
	ALL 				// 	3 : Print information about both virtual & physical Particles
};

enum class ParticleApproximation {
	NONE, 							// 	0	: No Particle Approximation
	Algorithm_1,					//	1 	: (e.g. ( p(i) + p(j) ) / ( rho(i) * rho(j) )
	Algorithm_2						//	2 	: (e.g. ( p(i) / rho(i)^2 + p(j) / rho(j)^2 )
};

enum class NearestNeighbourSearch {
	NONE, 							//	0	: No Particle Search Algorithm
	DirectSearch,					// 	1 	: Simplest and direct searching
	LinkedList,						// 	2 	: Sorting grid linked list
	TreeSearch						// 	3 	: Tree algorithm
};

enum class SmoothingLengthEvolution {
	NONE, 							//	0 	: Keep unchanged
	DensityAlgebraic,				//	1	:  							h 		= fac * ( m / rho ) * ( 1 / dim )
	SmoothingLengthODE,				//	2	: 		 					dh/dt	= ( -1 / dim ) * ( h / rho ) * ( drho / dt )
	OtherApproach 					//	3 	: Other approaches ( e.g.	h		= h_0 * ( rho_0 / rho )^( 1 / dim ) )
};

enum class SmoothingKernel {
	NONE, 							// 	0	: No Kernel Function
	CubicSpline,					// 	1	: Cubic Spline Kernel Function by W4 - Spline (Monaghan 1985)
	Gauss,							// 	2	: Gauss Kernel Function (Gingold and Monaghan 1981)
	Quintic 						// 	3	: Quintic Kernel Function (Morris 1997)
};

enum class DensityCalculation {
	NONE, 							// 	0	: No Formula chosen
	SummationDensity, 				// 	1 	: Use sph formula on density
	CSPM, 							// 	2 	: Summation Density with density normalization by using CSPM
	ContinuityEquation				//	3 	: Use continuity equation
};

enum class VelocityAveraging {
	NONE,							//	0	: No average treatment.
	Monaghan						//	1 	: Monaghan treatment on average velocity
};

// Use for physical & virtual Particles
enum class InitialConfiguration {
	NONE, 							// 	0 	: No initial Configuration
	Generate,						//	1	: Generate initial configuration
	Load							// 	2	: Load initial configuration data
};

enum class LoadInitial {
	NONE, 							// 	0 	: Dont load anything
	Physical,						//	1	: Load initial physical particles
	Virtual							// 	2	: Load initial virtual particles
};

enum class ParticleType {
	NONE,							//	0	: No use of virtual particle
	Physical,						// 	1 	: Physical Particle
	Virtual							// 	2 	: Virtual particle
};

enum class PDEs {
	NONE,							// 	0	: Don't solve any PDEs
	EulerEquations,					//	1 	: Solve Euler Equations
	NavierStokesEquations			// 	2	: Solve Navier-Stokes Equations
};

enum class ExternalForce {
	NONE,							//	0 	: No external force
	Gravity, 						// 	1 	: Consider Gravity
	BoundaryForce 					// 	2 	: Consider Boundary particle & penalty anti-penetration force
};

enum class Artificial {
	NONE, 							//	0 	: No artificial Quantity
	Viscosity,						// 	1 	: Consider artificial viscosity,
	Heat							// 	2 	: Consider artificial heating,
};

enum class Symmetry {
	NONE, 							//	0 	: no symmetry
	Axis, 							//	1 	: axis symmetry
	Center							// 	2 	: center symmetry
};

enum class ProblemType {
	NONE,
	ShockTube, 						//	1	: carry out shock tube simulation
	ShearCavity 					// 	2	: carry out shear cavity simulation
};


// Some operators

// std::vector

template <class T>
std::vector<T> operator+ ( const std::vector<T> & v1, const std::vector<T> & v2   ){

	size_t n1 = v1.size();
	size_t n2 = v2.size();

	MYASSERT (n1 == n2, "Error trying to add two std::vectors! Size not equal!");

	std::vector<T> sum( n1, .0);
	for( size_t i = 0; i < n1; i++){
		sum[i] =  v1[i] + v2[i];
	}

	return sum;
}

template <class T>
std::vector<T> & operator+= ( std::vector<T> & v1, const std::vector<T> & v2   ){

	return v1 = v1 + v2;
}

template <class T>
std::vector<T> operator* ( const std::vector<T> & v, const T & s   ){

	size_t n = v.size();

	std::vector<T> pro( n, .0);
	for( size_t i = 0; i < n; i++){
		pro[i] =  s*v[i];
	}

	return pro;
}

template <class T>
std::vector<T> operator* ( const T & s, const std::vector<T> & v  ){

	return v*s;
}

template <class T>
std::vector<T>& operator*= ( std::vector<T> & v, const T & s  ){

	return v = v*s;
}

} /* namespace SPH */

#endif /* PRIMITIVES_DATATYPES_MYTYPES_H_ */




























