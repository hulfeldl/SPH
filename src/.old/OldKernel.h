/*
 * Kernel.h
 *
 *  Created on: Aug 5, 2020
 *      Author: hulfeldl
 */




#ifndef KERNELS_OLDKERNEL_H_
#define KERNELS_OLDKERNEL_H_

#include "Primitives/DataTypes/MyTypes.h"
#include "Settings/Settings.h"
#include "Primitives/Constants.h"

#include <array>
#include <math.h>
#include <iostream>
#include <stdlib.h>

namespace SPH {

enum class kMember {
	DIM,
	kType,
	sle,
	hsml
};

template<uint8_t DIM , ParticleType pType >
class Particle;

template < uint8_t DIM , SmoothingKernel kType >
class AbstractKernel {

private:

	Particle<DIM> & m_Particle;


public:

	// Empty Constructor
	AbstractKernel();

	// Std. Constructor
	AbstractKernel( const int sle,
			const double hsml = 0 );

	// Std. Constructor
	AbstractKernel( const SmoothingLengthEvolution sle,
			const double hsml = 0 );

	virtual ~AbstractKernel();

	double hsml () const;

	void evolveSL();

	void update_hsml( const double hsml );

	void update_hsml(const double dhsmldt ,	const double dt );

};


template < uint8_t DIM , SmoothingKernel kType >
class OldKernel : public AbstractKernel<DIM,kType> {



public:

	// Empty Constructor
	OldKernel();

	// Std. Constructor
	OldKernel( const int sle,
			const double hsml = 0 );

	// Std. Constructor
	OldKernel( const SmoothingLengthEvolution sle,
			const double hsml = 0 );

	virtual ~OldKernel();

	constexpr double W() const;

	double W( const double r , const v4df & dx ) const;

	constexpr v4df dWdx() const;

	v4df dWdx( const double r , const v4df & dx ) const;

	template <typename T>
	T get( const kMember & Member ) const;

private:

};

template < uint8_t DIM >
class OldKernel<DIM,SmoothingKernel::CubicSpline> : public AbstractKernel<DIM,SmoothingKernel::CubicSpline> {

	SmoothingLengthEvolution m_sle;	// 	sle 	: Smoothing Length Evolution Algorithm
	double m_hsml;					//	hsml 	: smoothing lengths of particles


public:

	constexpr double W() const;

	double W( const double r ) const;

	constexpr v4df dWdx() const;

	v4df dWdx( const double r , const v4df & dx ) const;

private:

	double factor () const;

};

template < uint8_t DIM >
class OldKernel<DIM,SmoothingKernel::Gauss> {

	SmoothingLengthEvolution m_sle;	// 	sle 	: Smoothing Length Evolution Algorithm
	double m_hsml;					//	hsml 	: smoothing lengths of particles


public:

	double W() const;

	double W( const double r ) const;

	constexpr v4df dWdx() const;

	v4df dWdx( const double r , const v4df & dx ) const;

};

template < uint8_t DIM >
class OldKernel<DIM,SmoothingKernel::Quintic> {

private:

	SmoothingLengthEvolution m_sle;	// 	sle 	: Smoothing Length Evolution Algorithm
	double m_hsml;					//	hsml 	: smoothing lengths of particles



public:

	double W() const;

	double W( const double r ) const;

	constexpr v4df dWdx() const;

	v4df dWdx( const double r , const v4df & dx ) const;

private:

	double factor () const;

};

} /* namespace SPH */

#endif /* KERNELS_OLDKERNEL_H_ */










