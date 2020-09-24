/*
 * Symmetric3x3.cpp
 *
 *  Created on: Aug 11, 2020
 *      Author: hulfeldl
 */

#include "Symmetric3x3.h"

namespace SPH {

Symmetric3x3::Symmetric3x3() :
	m_xx ( .0 ),
	m_xy ( .0 ),
	m_xz ( .0 ),
	m_yy ( .0 ),
	m_yz ( .0 ),
	m_zz ( .0 )
	{

}

Symmetric3x3::Symmetric3x3( const v6df & Entries ) :
	m_xx ( Entries[0] ),
	m_xy ( Entries[1] ),
	m_xz ( Entries[2] ),
	m_yy ( Entries[3] ),
	m_yz ( Entries[4] ),
	m_zz ( Entries[5] ) {

}

Symmetric3x3::Symmetric3x3( const double xx,
							const double xy,
							const double xz,
							const double yy,
							const double yz,
							const double zz ) :
	m_xx ( xx ),
	m_xy ( xy ),
	m_xz ( xz ),
	m_yy ( yy ),
	m_yz ( yz ),
	m_zz ( zz ) {

}

Symmetric3x3::~Symmetric3x3() = default;

double Symmetric3x3::xx() const {
	return m_xx;
}

double Symmetric3x3::xy() const {
	return m_xy;
}

double Symmetric3x3::xz() const {
	return m_xz;
}

double Symmetric3x3::yy() const {
	return m_yy;
}

double Symmetric3x3::yz() const {
	return m_yz;
}

double Symmetric3x3::zz() const {
	return m_zz;
}

Symmetric3x3 Symmetric3x3::operator+	( const Symmetric3x3 & s ) const {

	Symmetric3x3 sum;

	sum.m_xx 	= this->m_xx + s.m_xx;
	sum.m_xx 	= this->m_xx + s.m_xx;
	sum.m_xx 	= this->m_xx + s.m_xx;
	sum.m_xx 	= this->m_xx + s.m_xx;
	sum.m_xx 	= this->m_xx + s.m_xx;
	sum.m_xx 	= this->m_xx + s.m_xx;

	return sum;
}

Symmetric3x3& Symmetric3x3::operator+=	( const Symmetric3x3 & s ) {

	*this	= *this + s;

//	this->m_xx += s.m_xx;
//	this->m_xy += s.m_xy;
//	this->m_xz += s.m_xz;
//	this->m_yy += s.m_yy;
//	this->m_yz += s.m_yz;
//	this->m_zz += s.m_zz;

	return *this;
}

Symmetric3x3& Symmetric3x3::operator= ( const Symmetric3x3 & s ) {

    // self-assignment guard
    if ( this == &s )
        return *this;

    // do the copy
    this->m_xx = s.m_xx;
    this->m_xy = s.m_xy;
    this->m_xz = s.m_xz;
    this->m_yy = s.m_yy;
    this->m_yz = s.m_yz;
    this->m_zz = s.m_zz;

    // return the existing object so we can chain this operator
    return *this;
}

} /* namespace SPH */





















