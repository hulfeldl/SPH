/*
 * MyVector.h
 *
 *  Created on: Aug 13, 2020
 *      Author: hulfeldl
 */




#ifndef PRIMITIVES_DATATYPES_MYVECTOR_H_
#define PRIMITIVES_DATATYPES_MYVECTOR_H_

#include <array>
#include <vector>
#include <iostream>
#include <initializer_list>

#include "Primitives/MyError.h"

namespace SPH {

template < typename T , uint8_t DIM >
class Symmetric;

template < typename T , uint8_t DIM  >
class MyVector {

	typedef std::array < T , DIM > 	vNT;

private:

	vNT m_Data;

public:

	MyVector();
	template < typename... Args >
	MyVector( Args... args ) {

		const size_t nargs = sizeof...(Args);

//	    std::cerr << "n_args: " << std::to_string(nargs) << ", DIM: " << std::to_string(DIM) << std::endl;
//	    ( std::cout << ... << args ) << " " << std::endl;

	    MYASSERT( nargs == DIM, "Wrong number of input arguments!");

	    m_Data = std::array<std::common_type_t<Args...>, DIM>{args...};

	}

	MyVector( const MyVector & v );
	MyVector( const std::vector<T> & v );
	MyVector( const std::array<T,DIM> & v );

	virtual ~MyVector();

	T& operator[] ( uint8_t i);

	const T& operator[] ( uint8_t i) const;

	T operator^ ( const MyVector & v ) const; 	// Dot Product

	MyVector operator+ ( const MyVector & v ) const;

	MyVector operator- () const;

	MyVector operator- ( const MyVector & v ) const;

	MyVector& operator+= ( const MyVector & v );

	MyVector& operator-= ( const MyVector & v );

	MyVector& operator*= ( const T & s );

	template < typename TT , uint8_t D  >
	friend MyVector<TT,D> operator+ ( const MyVector<TT,D> & v1, const std::vector<TT> & v2 );

	template < typename TT , uint8_t D  >
	friend MyVector<TT,D> operator+ ( const std::vector<TT> & v1, const MyVector<TT,D> & v2 );

	template < typename TT , uint8_t D  >
	friend MyVector<TT,D> operator* ( const MyVector<TT,D> & v1, const std::vector<TT> & v2 );

	template < typename TT , uint8_t D  >
	friend MyVector<TT,D> operator* ( const std::vector<TT> & v1, const MyVector<TT,D> & v2 );

	template < typename TT , uint8_t D  >
	friend MyVector<TT,D> operator* ( const MyVector<TT,D> & v, const TT & s );

	template < typename TT , uint8_t D  >
	friend MyVector<TT,D> operator* ( const TT & s, const MyVector<TT,D> & v );

	template < typename TT , uint8_t D  >
	friend Symmetric<TT,D> h ( const MyVector<TT,D> & v1, const MyVector<TT,D> & v2 );

	template < typename TT , uint8_t D  >
	friend  std::ostream & operator<< ( std::ostream & s, const MyVector<TT,D> & v );
};

} /* namespace SPH */

#endif /* PRIMITIVES_DATATYPES_MYVECTOR_H_ */


















