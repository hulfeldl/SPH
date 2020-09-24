/*
 * MyMatrix.cpp
 *
 *  Created on: Aug 13, 2020
 *      Author: hulfeldl
 */




#include <cstdarg>
#include <type_traits>

#include "MyMatrix.h"

namespace SPH {

template < typename T , uint8_t DIM1 , uint8_t DIM2  >
MyMatrix<T,DIM1,DIM2>::MyMatrix() {
	// TODO Auto-generated constructor stub

}

template < typename T , uint8_t DIM1 , uint8_t DIM2  >
MyMatrix<T,DIM1,DIM2>::MyMatrix( const vNT & n_args... ) {

    std::va_list args;
    va_start(args, n_args);

    MYASSERT(n_args == DIM1, "Wrong number of input arguments!");

    for (size_t i = 0; i < n_args; i++ ){
    	m_Data[i] = va_arg( args,vNT);
    }

    va_end(args);

}

template < typename T , uint8_t DIM1 , uint8_t DIM2  >
MyMatrix<T,DIM1,DIM2>::~MyMatrix() {
	// TODO Auto-generated destructor stub
}

template < typename T , uint8_t DIM1, uint8_t DIM2 >
MyMatrix<T,DIM1,DIM2>& MyMatrix<T,DIM1,DIM2>::operator= ( const MyMatrix<T,DIM1,DIM2> & m ) {

    // self-assignment guard
    if ( this == &m )
        return *this;

    m_Data = m.m_Data;

    // return the existing object so we can chain this operator
    return *this;
}

template < typename T , uint8_t DIM1 , uint8_t DIM2  >
MyVector < T , DIM2 > & MyMatrix<T,DIM1,DIM2>::operator[] ( uint8_t i){
	return m_Data[i];
}

template < typename T , uint8_t DIM1, uint8_t DIM2 >
MyMatrix<T,DIM1,DIM2> MyMatrix<T,DIM1,DIM2>::operator+	( const MyMatrix<T,DIM1,DIM2> & m ) const {

	MyMatrix<T,DIM1,DIM2> sum;
	for ( uint8_t i = 0; i < DIM1; i++ ){
		sum[i] = m_Data[i] + m.m_Data[i];
	}

	return sum;
}

template < typename T , uint8_t DIM1, uint8_t DIM2 >
MyMatrix<T,DIM1,DIM2>& MyMatrix<T,DIM1,DIM2>::operator+= ( const MyMatrix<T,DIM1,DIM2> & m ) {

	return *this = *this + m;
}

template < typename T , uint8_t DIM1 , uint8_t DIM2  >
T MyMatrix<T,DIM1,DIM2>::operator^ ( const MyMatrix<T,DIM1,DIM2> & m ) const {

	T sum = T();
	for ( uint8_t i = 0; i < DIM1; i++ ){
		sum += m_Data[i] ^ m_Data[i];
	}

	return sum;

}

//template class MyMatrix<double,1,1>;

} /* namespace SPH */



















