/*
 * QuinticKernel.h
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */

#ifndef KERNELS_QUINTICKERNEL_H_
#define KERNELS_QUINTICKERNEL_H_

#include "Kernel.h"

namespace SPH {

template < uint8_t DIM >
class QuinticKernel: public Kernel<DIM> {

public:

	QuinticKernel() : Kernel<DIM>(){}
	QuinticKernel( uint8_t sle, double hsml ) : Kernel<DIM>(sle,hsml){}
	QuinticKernel( SmoothingLengthEvolution sle, double hsml ) : Kernel<DIM>(sle,hsml){}
	virtual ~QuinticKernel();
	QuinticKernel(const QuinticKernel &other);
	QuinticKernel(QuinticKernel &&other);
	QuinticKernel& operator=(const QuinticKernel &other);
	QuinticKernel& operator=(QuinticKernel &&other);

	constexpr double k() const;

	double W() const;

	double W( const double r ) const;

	constexpr vDdf<DIM> dWdx() const;

	vDdf<DIM> dWdx( const double r , const vDdf<DIM> & dx ) const;

private:

	double factor () const;

};

} /* namespace SPH */

#endif /* KERNELS_QUINTICKERNEL_H_ */
