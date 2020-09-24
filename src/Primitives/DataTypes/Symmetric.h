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

#include "MyMatrix.h"

namespace SPH {

template < typename T , uint8_t DIM >
class Symmetric  {

	enum k : uint32_t {
		Diag 	= DIM,				// Number of Diagonal Elements
		Other	= DIM*(DIM-1)/2,	// Number of other unique Elements
		Unique 	= DIM*(DIM+1)/2		// Number of unique Elements
	};

	typedef MyVector<T,k::Diag> 		vDiagT;
	typedef MyVector<T,k::Other> 		vOthT;

private:

	vDiagT 	m_Diag;
	vOthT 	m_Other;

public:

	Symmetric();

	Symmetric(	const vDiagT & Diag, const vOthT & Other );

	Symmetric( 	const std::initializer_list<T> &  Entries );

	Symmetric( 	const T Entries... );

	virtual ~Symmetric();

	Symmetric& 	operator=	( const Symmetric & sym );

	MyVector<T,DIM> operator[] ( const uint8_t i) const;

	T & operator() ( const uint8_t i, const uint8_t j );

	const T & operator() ( const uint8_t i, const uint8_t j ) const;

	Symmetric 	operator+	( const Symmetric & sym ) const;

	MyVector<T,DIM> operator* ( const MyVector<T,DIM> & v ) const;

	Symmetric& 	operator+=	( const Symmetric & sym );

	T operator^ ( const Symmetric & sym ) const; 	// Element wise multiplication

	template < typename TT , uint8_t D >
	friend Symmetric<TT,D> operator* ( const Symmetric<TT,D>  & sym, const TT & s );

	template < typename TT , uint8_t D >
	friend Symmetric<TT,D>  operator* ( const TT & s, const Symmetric<TT,D>  & sym );


private:

	static uint32_t index( const uint8_t i, const uint8_t j);
};



// Partial specializations
// ------------------------------------------------------------------------------------------
//template < typename T >
//class Symmetric<T,1> : public MyMatrix<T,1,1> {
//
//private:
//
//public:
//
//	Symmetric();
//
//	Symmetric(	const Symmetric<T,1> & s );
//
//	Symmetric( 	const T xx );
//
//	virtual ~Symmetric();
//
//	T xx() const;
//
//};
//
//template < typename T >
//class Symmetric<T,2> : public MyMatrix<T,2,2> {
//
//private:
//
//public:
//
//	Symmetric();
//
//	Symmetric(	const Symmetric<T,2> & s );
//
//	Symmetric( 	const T xx,
//				const T xy,
//				const T yy );
//
//	virtual ~Symmetric();
//
//	T xx() const;
//	T xy() const;
//	T yy() const;
//
//};
//
//template < typename T >
//class Symmetric<T,3> : public MyMatrix<T,3,3> {
//
//private:
//
//public:
//
//	Symmetric();
//
//	Symmetric(	const Symmetric<T,3> & s );
//
//	Symmetric( 	const T xx,
//				const T xy,
//				const T xz,
//				const T yy,
//				const T yz,
//				const T zz );
//
//	virtual ~Symmetric();
//
//	T xx() const;
//	T xy() const;
//	T xz() const;
//	T yy() const;
//	T yz() const;
//	T zz() const;
//
//};

} /* namespace SPH */

#endif /* PRIMITIVES_SYMMETRIC3X3_H_ */













