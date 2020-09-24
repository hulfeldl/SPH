/*
 * Symmetric3x3.h
 *
 *  Created on: Aug 11, 2020
 *      Author: hulfeldl
 */




#ifndef PRIMITIVES_SYMMETRIC3X3_H_
#define PRIMITIVES_SYMMETRIC3X3_H_

#include <array>
#include <vector>

namespace SPH {

class Symmetric3x3 {

	// Typedefs
	typedef std::array<double,4> 	v4df;
	typedef std::array<double,6> 	v6df;
	typedef std::vector<double> 	vNdf;
	typedef std::vector<v4df> 		vN4df;

private:

	double m_xx,	m_xy, 	m_xz,
					m_yy,	m_yz,
							m_zz;

public:

	Symmetric3x3();

	Symmetric3x3( const v6df & Entries );

	Symmetric3x3( 	const double xx,
					const double xy,
					const double xz,
					const double yy,
					const double yz,
					const double zz );

	virtual ~Symmetric3x3();

	double xx() const;
	double xy() const;
	double xz() const;
	double yy() const;
	double yz() const;
	double zz() const;

	Symmetric3x3 	operator+	( const Symmetric3x3 & s ) const;

	Symmetric3x3& 	operator+=	( const Symmetric3x3 & s );

	Symmetric3x3& 	operator=	( const Symmetric3x3 & s );

};

} /* namespace SPH */

#endif /* PRIMITIVES_SYMMETRIC3X3_H_ */













