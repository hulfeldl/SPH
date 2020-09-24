/*
 * MyVector.cpp
 *
 *  Created on: Aug 13, 2020
 *      Author: hulfeldl
 */




#include <cstdarg>

#include "MyVector.h"
#include "Symmetric.h"

namespace SPH {

template < typename T , uint8_t DIM  >
MyVector<T,DIM>::MyVector() :
	m_Data() {

}

template < typename T , uint8_t DIM  >
MyVector<T,DIM>::MyVector( const MyVector<T,DIM> & v ) :
	m_Data(v.m_Data) {

}

template < typename T , uint8_t DIM  >
MyVector<T,DIM>::MyVector( const std::vector<T> & v ) {

	MYASSERT( 	DIM == v.size(),
				"Failed to initialize MyVector with std::vector! Dimensions do not agree!" );

	for( uint8_t i = 0; i < DIM; i++ ) {
		m_Data[i] = v[i];
	}

}

template < typename T , uint8_t DIM  >
MyVector<T,DIM>::MyVector( const std::array<T,DIM> & v ) :
	m_Data(v) {
}


template < typename T , uint8_t DIM  >
MyVector<T,DIM>::~MyVector() = default;

template < typename T , uint8_t DIM  >
T& MyVector<T,DIM>::operator[] ( uint8_t i) {
	return m_Data[i];
}

template < typename T , uint8_t DIM  >
const T& MyVector<T,DIM>::operator[] ( uint8_t i) const {
	return m_Data[i];
}

template < typename T , uint8_t DIM  >
T MyVector<T,DIM>::operator^ ( const MyVector<T,DIM> & v ) const {

	T sum = T();
	for ( uint8_t i = 0; i < DIM; i++ ){
		sum += m_Data[i] * v.m_Data[i];
	}

	return sum;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> MyVector<T,DIM>::operator+ ( const MyVector<T,DIM> & v ) const {

	MyVector<T,DIM> sum;
	for ( uint8_t i = 0; i < DIM; i++ ){
		sum[i] = m_Data[i] + v.m_Data[i];
	}

	return sum;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> MyVector<T,DIM>::operator- () const {

	MyVector<T,DIM> diff;
	for ( uint8_t i = 0; i < DIM; i++ ){
		diff[i] = -m_Data[i];
	}

	return diff;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> MyVector<T,DIM>::operator- ( const MyVector<T,DIM> & v ) const {

	MyVector<T,DIM> diff;
	for ( uint8_t i = 0; i < DIM; i++ ){
		diff[i] = m_Data[i] - v.m_Data[i];
	}

	return diff;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> & MyVector<T,DIM>::operator+= ( const MyVector<T,DIM> & v ) {
	return *this = *this + v;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> & MyVector<T,DIM>::operator-= ( const MyVector<T,DIM> & v ) {
	return *this = *this - v;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> & MyVector<T,DIM>::operator*= ( const T & s ) {
	return *this = *this * s;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> operator+ ( const MyVector<T,DIM> & v1, const std::vector<T> & v2 ) {

	MYASSERT( 	DIM == v2.size(),
				"Failed to add MyVector with std::vector! Dimensions do not agree!" );

	MyVector<T,DIM> sum;
	for ( uint8_t i = 0; i < DIM; i++ ){
		sum[i] = v1.m_Data[i] + v2[i];
	}

	return sum;
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> operator+ ( const std::vector<T> & v1, const MyVector<T,DIM> & v2 ) {
	return operator+(v2, v1);
}

template < typename T , uint8_t DIM  >
MyVector<T,DIM> operator* ( const MyVector<T,DIM> & v, const T & s ) {

	MyVector<T,DIM> prod;

	for ( uint8_t i = 0; i < DIM; i++ ){
		prod[i] = s * v.m_Data[i];
	}

	return prod;

}

template < typename T, uint8_t DIM  >
MyVector<T,DIM> operator* ( const T & s, const MyVector<T,DIM> & v ) {
	return v * s;
}

template < typename T , uint8_t D  >
Symmetric<T,D> h ( const MyVector<T,D> & v1, const MyVector<T,D> & v2 ) {

	MyVector<T,D> 			Diag;
	MyVector<T,D*(D-1)/2> 	Other;

	uint32_t iO = 0;
	for ( uint8_t i = 0; i < D; i++ ){
		for ( uint8_t j = 0; i < D; i++ ){
			if ( i == j ){
				Diag[i] += static_cast<T>(2.0) * v1[j] * v2[j];
			}
			else {
				Diag[i] -= v1[j] * v2[j];
			}
		}

    	for ( uint8_t j = i+1; j < D; j++ ){
    		Other[iO] = v1[i] * v2[j] + v1[j] * v2[i];
    		iO++;
    	}
	}
	Diag *= static_cast<T>(2.0) / static_cast<T>(3.0);

	return Symmetric<T,D>( Diag, Other );
}


template < typename TT , uint8_t DIM  >
std::ostream & operator<< ( std::ostream & os, const MyVector<TT,DIM> & v ) {

	for (int d = 0; d < DIM; d++ ) {
		os << v[d] << " ";

	}

	return os;
}


template class MyVector<double,0>;
template class MyVector<double,1>;
template class MyVector<double,2>;
template class MyVector<double,3>;


template MyVector<double,1> operator+ ( const MyVector<double,1> & v1, const std::vector<double> & v2 );
template MyVector<double,2> operator+ ( const MyVector<double,2> & v1, const std::vector<double> & v2 );
template MyVector<double,3> operator+ ( const MyVector<double,3> & v1, const std::vector<double> & v2 );

template MyVector<double,1> operator+ ( const std::vector<double> & v1, const MyVector<double,1> & v2 );
template MyVector<double,2> operator+ ( const std::vector<double> & v1, const MyVector<double,2> & v2 );
template MyVector<double,3> operator+ ( const std::vector<double> & v1, const MyVector<double,3> & v2 );

//template MyVector<double,1> operator* ( const MyVector<double,1> & v1, const std::vector<double> & v2 );
//template MyVector<double,2> operator* ( const MyVector<double,2> & v1, const std::vector<double> & v2 );
//template MyVector<double,3> operator* ( const MyVector<double,3> & v1, const std::vector<double> & v2 );
//
//template MyVector<double,1> operator* ( const std::vector<double> & v1, const MyVector<double,1> & v2 );
//template MyVector<double,2> operator* ( const std::vector<double> & v1, const MyVector<double,2> & v2 );
//template MyVector<double,3> operator* ( const std::vector<double> & v1, const MyVector<double,3> & v2 );

template MyVector<double,0> operator* ( const MyVector<double,0> & v, const double & s );
template MyVector<double,1> operator* ( const MyVector<double,1> & v, const double & s );
template MyVector<double,2> operator* ( const MyVector<double,2> & v, const double & s );
template MyVector<double,3> operator* ( const MyVector<double,3> & v, const double & s );

template MyVector<double,0> operator* ( const double & s, const MyVector<double,0> & v );
template MyVector<double,1> operator* ( const double & s, const MyVector<double,1> & v );
template MyVector<double,2> operator* ( const double & s, const MyVector<double,2> & v );
template MyVector<double,3> operator* ( const double & s, const MyVector<double,3> & v );

template Symmetric<double,1> h ( const MyVector<double,1> & v1, const MyVector<double,1> & v2 );
template Symmetric<double,2> h ( const MyVector<double,2> & v1, const MyVector<double,2> & v2 );
template Symmetric<double,3> h ( const MyVector<double,3> & v1, const MyVector<double,3> & v2 );

template std::ostream & operator<< ( std::ostream & os, const MyVector<double,1> & v );
template std::ostream & operator<< ( std::ostream & os, const MyVector<double,2> & v );
template std::ostream & operator<< ( std::ostream & os, const MyVector<double,3> & v );

} /* namespace SPH */














































