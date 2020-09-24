/*
 * Kernel.h
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */




#ifndef KERNELS_KERNEL_H_
#define KERNELS_KERNEL_H_

#include <math.h>

#include "Primitives/DataTypes/MyTypes.h"
#include "Settings/Settings.h"
#include "Primitives/Constants.h"

namespace SPH {

constexpr uint8_t ALIGN_LENGTH = 32;

template < uint8_t DIM >
class alignas(ALIGN_LENGTH) Kernel {

	typedef std::unique_ptr<Kernel> 	kPtr;

private:

	SmoothingLengthEvolution m_sle;	// 	sle 	: Smoothing Length Evolution Algorithm
	double m_hsml;					//	hsml 	: smoothing lengths of particles
	double m_f; 					//	f 		: Factor for calculating W and dWdx

public:

	Kernel() : m_sle(), m_hsml(), m_f( ){}
	Kernel(uint8_t sle, double hsml) : m_sle(static_cast<SmoothingLengthEvolution>(sle)), m_hsml(hsml), m_f(){ }
	Kernel(SmoothingLengthEvolution sle, double hsml) : m_sle(sle), m_hsml(hsml), m_f( ){}
	virtual ~Kernel() = default;
	Kernel(const Kernel &other);
	Kernel(Kernel &&other);
	Kernel& operator=(const Kernel &other);
	Kernel& operator=(Kernel &&other);

	virtual kPtr clone() const = 0;

	double hsml () const { return m_hsml; }

	void update_hsml( const double hsml ){ m_hsml = hsml; }

	void update_hsml( const double dhdt, const double dt ){ m_hsml += dhdt*dt; }

	SmoothingLengthEvolution sle() const { return m_sle; }

	virtual constexpr double k() const = 0;

	virtual void update_f() = 0;

	virtual double W() const = 0;

	virtual double W( const double r ) const = 0;

	virtual constexpr vDdf<DIM> dWdx() const = 0;

	virtual vDdf<DIM> dWdx( const double r , const vDdf<DIM> & dx ) const = 0;

protected:

	double f() const { return m_f; }

	void set_f( const double f ) { m_f = f; }

private:

	virtual double factor () const = 0;

};

} /* namespace SPH */

#endif /* KERNELS_KERNEL_H_ */















