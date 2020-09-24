/*
 * GaussKernel.h
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */

#ifndef KERNELS_GAUSSKERNEL_H_
#define KERNELS_GAUSSKERNEL_H_

#include "Kernel.h"

namespace SPH {

template < uint8_t DIM >
class GaussKernel: public Kernel<DIM> {

public:

	GaussKernel() : Kernel<DIM>(){}
	GaussKernel( uint8_t sle, double hsml ) : Kernel<DIM>(sle,hsml){}
	GaussKernel( SmoothingLengthEvolution sle, double hsml ) : Kernel<DIM>(sle,hsml){}
	virtual ~GaussKernel();
	GaussKernel(const GaussKernel &other);
	GaussKernel(GaussKernel &&other);
	GaussKernel& operator=(const GaussKernel &other);
	GaussKernel& operator=(GaussKernel &&other);

	constexpr double k() const;

	double W() const;

	double W( const double r ) const;

	constexpr vDdf<DIM> dWdx() const;

	vDdf<DIM> dWdx( const double r , const vDdf<DIM> & dx ) const;

private:

	double factor () const;

};

} /* namespace SPH */

#endif /* KERNELS_GAUSSKERNEL_H_ */
