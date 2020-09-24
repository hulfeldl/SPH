/*
 * Symmetric3x3.cpp
 *
 *  Created on: Aug 11, 2020
 *      Author: hulfeldl
 */




#include <cstdarg>
#include <algorithm>

#include "Symmetric.h"

namespace SPH {

template < typename T , uint8_t DIM >
Symmetric<T,DIM>::Symmetric() = default;

template < typename T , uint8_t DIM >
Symmetric<T,DIM>::Symmetric(const vDiagT & Diag, const vOthT & Other ) :
	m_Diag(Diag),
	m_Other(Other) {

}

template < typename T , uint8_t DIM >
Symmetric<T,DIM>::Symmetric( const std::initializer_list<T> &  Entries ) {

}

template < typename T , uint8_t DIM >
Symmetric<T,DIM>::Symmetric( const T n_args... ) {

    std::va_list args;
    va_start(args, n_args);

    MYASSERT(n_args == k::Unique, "Wrong number of input arguments!");

    uint32_t iO = 0;
    uint8_t nOther = DIM;
    for (uint8_t i = 0; i < DIM; i++ ){
    	m_Diag[i] = va_arg( args, T);
    	nOther--;
    	for ( uint8_t j = 0; j < nOther; j++ ){
    		m_Other[iO] = va_arg( args, T);
    		iO++;
    	}
    }

    va_end(args);
}


template < typename T , uint8_t DIM >
Symmetric<T,DIM>::~Symmetric() = default;

template < typename T , uint8_t DIM >
Symmetric<T,DIM>& 	Symmetric<T,DIM>::operator=	( const Symmetric<T,DIM> & sym ) {
	m_Diag 	= sym.m_Diag;
	m_Other = sym.m_Other;

	return *this;
}

template < typename T , uint8_t DIM >
MyVector<T,DIM> Symmetric<T,DIM>::operator[] ( const uint8_t i) const {

	MyVector<T,DIM> row;

	for (uint8_t j = 0; j < DIM; j++ ){
		row[j] = operator()(i,j);
	}

	return row;
}

template < typename T , uint8_t DIM >
T & Symmetric<T,DIM>::operator() ( const uint8_t i, const uint8_t j ) {

	if ( i == j ) {
		return m_Diag[i];
	}
	else {
		return m_Other[index(i,j)];
	}
}

template < typename T , uint8_t DIM >
const T & Symmetric<T,DIM>::operator() ( const uint8_t i, const uint8_t j ) const {

	if ( i == j ) {
		return m_Diag[i];
	}
	else {
		return m_Other[index(i,j)];
	}
}

template < typename T , uint8_t DIM >
uint32_t Symmetric<T,DIM>::index( const uint8_t i, const uint8_t j) {

	if ( i < j ){
		return i*(DIM* - (i+1)/2) + (j - (i+1) );
	}
	else {
		return j*(DIM* - (j+1)/2) + (i - (j+1) );
	}

}

template < typename T , uint8_t DIM >
Symmetric<T,DIM> Symmetric<T,DIM>::operator+ ( const Symmetric<T,DIM> & sym ) const {

	Symmetric<T,DIM> sum;

	sum.m_Diag 	= m_Diag + sym.m_Diag;
	sum.m_Other = m_Other + sym.m_Other;

	return sum;
}

template < typename T , uint8_t DIM >
MyVector<T,DIM> Symmetric<T,DIM>::operator* ( const MyVector<T,DIM> & v ) const {

	MyVector<T,DIM> pro;

	for( uint8_t i = 0; i < DIM; i++ ){
		pro[i] = operator[](i) ^ v;
	}

	return pro;
}

template < typename T , uint8_t DIM >
Symmetric<T,DIM>& 	Symmetric<T,DIM>::operator+= ( const Symmetric<T,DIM> & sym ) {
	return *this = *this + sym;
}

// Element wise multiplication
template < typename T , uint8_t DIM >
T Symmetric<T,DIM>::operator^ ( const Symmetric<T,DIM> & sym ) const {
	T sumO = m_Other ^ m_Other;
	return (m_Diag ^ m_Diag) +  sumO + sumO;
}

template < typename T , uint8_t DIM >
Symmetric<T,DIM> operator* ( const Symmetric<T,DIM> & sym, const T & s ) {

	Symmetric<T,DIM> pro;

	pro.m_Diag 	= s * sym.m_Diag;
	pro.m_Other = s * sym.m_Other;

	return pro;
}

template < typename T , uint8_t DIM >
Symmetric<T,DIM> operator* ( const T & s, const Symmetric<T,DIM> & sym ) {
	return sym*s;
}


//
//
//// Partially Specialized 1D
////-------------------------------------------------------------------------------------------
//template < typename T >
//Symmetric<T,1>::Symmetric() {
//
//}
//
//template < typename T >
//Symmetric<T,1>::Symmetric(	const Symmetric<T,1> & s ){
//
//}
//
//template < typename T >
//Symmetric<T,1>::Symmetric( 	const T xx ) :
//	MyMatrix<T,1,1>( MyVector<T,1>( xx ) ) {
//
//}
//
//
//template < typename T >
//Symmetric<T,1>::~Symmetric() = default;
//
//template < typename T >
//T Symmetric<T,1>::xx() const {
//	return MyMatrix<T,1,1>::operator[](0)[0];
//}
//
//
//// Partially Specialized 2D
////-------------------------------------------------------------------------------------------
//
//template < typename T >
//Symmetric<T,2>::Symmetric();
//
//template < typename T >
//Symmetric<T,2>::Symmetric(	const Symmetric<T,2> & s );
//
//template < typename T >
//Symmetric<T,2>::Symmetric( 	const T xx,
//							const T xy,
//							const T yy ) :
//	m_Data( MyVector<T,2>( xx, xy ),
//			MyVector<T,2>( xy, yy ) ) {
//
//}
//
//template < typename T >
//Symmetric<T,2>::~Symmetric();
//
//template < typename T >
//T Symmetric<T,2>::xx() const {
//	return m_Data[0][0];
//}
//
//template < typename T >
//T Symmetric<T,2>::xy() const {
//	return m_Data[0][1];
//}
//
//template < typename T >
//T Symmetric<T,2>::yy() const {
//	return m_Data[1][1];
//}
//
//// Partially Specialized 3D
////-------------------------------------------------------------------------------------------
//
//template < typename T >
//Symmetric<T,3>::Symmetric();
//
//template < typename T >
//Symmetric<T,3>::Symmetric(	const Symmetric<T,3> & s );
//
//template < typename T >
//Symmetric<T,3>::Symmetric( 	const T xx,
//							const T xy,
//							const T xz,
//							const T yy,
//							const T yz,
//							const T zz ) :
//	m_Data( MyVector<T,3>( xx, xy, xz ),
//			MyVector<T,3>( xy, yy, yz ),
//			MyVector<T,3>( xz, yz, zz ) ) {
//
//}
//
//template < typename T >
//Symmetric<T,3>::~Symmetric();
//
//template < typename T >
//T Symmetric<T,3>::xx() const {
//	return m_Data[0][0];
//}
//
//template < typename T >
//T Symmetric<T,3>::xy() const {
//	return m_Data[0][1];
//}
//
//template < typename T >
//T Symmetric<T,3>::xz() const {
//	return m_Data[0][2];
//}
//
//template < typename T >
//T Symmetric<T,3>::yy() const {
//	return m_Data[1][1];
//}
//
//template < typename T >
//T Symmetric<T,3>::yz() const {
//	return m_Data[1][2];
//}
//
//template < typename T >
//T Symmetric<T,3>::zz() const {
//	return m_Data[2][2];
//}


template class Symmetric<double,1>;
template class Symmetric<double,2>;
template class Symmetric<double,3>;


template Symmetric<double,1> operator* ( const Symmetric<double,1>  & sym, const double & s );
template Symmetric<double,2> operator* ( const Symmetric<double,2>  & sym, const double & s );
template Symmetric<double,3> operator* ( const Symmetric<double,3>  & sym, const double & s );

template Symmetric<double,1>  operator* ( const double & s, const Symmetric<double,1>  & sym );
template Symmetric<double,2>  operator* ( const double & s, const Symmetric<double,2>  & sym );
template Symmetric<double,3>  operator* ( const double & s, const Symmetric<double,3>  & sym );

} /* namespace SPH */





















