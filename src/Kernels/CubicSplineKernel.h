/*
 * CubicSplineKernel.h
 *
 *  Created on: Sep 1, 2020
 *      Author: hulfeldl
 */

#ifndef KERNELS_CUBICSPLINEKERNEL_H_
#define KERNELS_CUBICSPLINEKERNEL_H_

#include "Kernel.h"

namespace SPH {

template < uint8_t DIM >
class CubicSplineKernel: public Kernel<DIM> {

	typedef std::unique_ptr<Kernel<DIM>> 	kPtr;

public:

	CubicSplineKernel() : Kernel<DIM>(){ update_f(); }
	CubicSplineKernel( uint8_t sle, double hsml ) : Kernel<DIM>(sle,hsml){ update_f(); }
	CubicSplineKernel( SmoothingLengthEvolution sle, double hsml ) : Kernel<DIM>(sle,hsml){ update_f();}
	virtual ~CubicSplineKernel();
	CubicSplineKernel(const CubicSplineKernel &other) : Kernel<DIM>( other.sle(), other.hsml() ){ update_f(); }
	CubicSplineKernel(CubicSplineKernel &&other);
	CubicSplineKernel& operator=(const CubicSplineKernel &other);
	CubicSplineKernel& operator=(CubicSplineKernel &&other);

	kPtr clone() const { return std::make_unique<CubicSplineKernel>(*this); }

	constexpr double k() const { return 2.0; }

	void update_f();

	double W() const;

	double W( const double r ) const;

	constexpr vDdf<DIM> dWdx() const { return vDdf<DIM>(); }

	vDdf<DIM> dWdx( const double r , const vDdf<DIM> & dx ) const;

private:

	double factor () const;

};

} /* namespace SPH */

#endif /* KERNELS_CUBICSPLINEKERNEL_H_ */
