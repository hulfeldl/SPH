/*
 * MyMatrix.h
 *
 *  Created on: Aug 13, 2020
 *      Author: hulfeldl
 */




#ifndef PRIMITIVES_MYMATRIX_H_
#define PRIMITIVES_MYMATRIX_H_

#include "MyVector.h"

namespace SPH {

template < typename T , uint8_t DIM1 , uint8_t DIM2  >
class MyMatrix {

	typedef MyVector < T , DIM2 > 		vNT;
	typedef std::array < vNT , DIM1 > 	mMNT;

private:

	mMNT m_Data;

public:

	MyMatrix();

	MyMatrix( const vNT & v... );

	virtual ~MyMatrix();

	MyMatrix& 	operator=	( const MyMatrix & m );

	vNT& operator[] ( uint8_t i);

	MyMatrix 	operator+	( const MyMatrix & m ) const;

	MyMatrix& 	operator+=	( const MyMatrix & m );

	T operator^ ( const MyMatrix & m ) const; 	// Element wise multiplication
};

} /* namespace SPH */

#endif /* PRIMITIVES_MYMATRIX_H_ */













