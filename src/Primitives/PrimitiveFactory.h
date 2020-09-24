/*
 * PrimitiveFactory.h
 *
 *  Created on: Aug 18, 2020
 *      Author: hulfeldl
 */




#ifndef PRIMITIVES_PRIMITIVEFACTORY_H_
#define PRIMITIVES_PRIMITIVEFACTORY_H_

#include <memory>
#include <map>
#include <vector>

#include "Primitives/MyError.h"

//template < typename T, class Base >
//class PrimitiveFactory {
//
//	typedef std::unique_ptr<Base>	BasePtr;
//	typedef std::map<T,BasePtr> 	mapType;
//
//
//private:
//
//	mapType m_FacMap;
//
//public:
//
//	template < class Derived >
//	PrimitiveFactory( const std::pair<T,Derived> & Args...  );
//
//	template < class Derived >
//	PrimitiveFactory( const T keys, const Derived & vals... );
//
//	virtual ~PrimitiveFactory();
//
//	BasePtr create ( const T Key );
//
//
//};

template < class Base , template<uint8_t> class Derived >
class DimFactory {

	typedef std::unique_ptr<Base>	BasePtr;

public:

	template< class... Args>
	BasePtr create ( uint8_t DIM, Args... args ) const {

		switch (DIM) {
		case 1:
			return std::make_unique<Derived<1>>(args...);

		case 2:
			return std::make_unique<Derived<2>>(args...);

		case 3:
			return std::make_unique<Derived<3>>(args...);

		default:
			MYASSERT (false, "Dimension must be 1, 2 or 3!");
		}

	}
};




// My Factory
// ------------------------------------------------------------------
//
//template < class Base >
//class Creator {
//
//	typedef std::unique_ptr<Base> basePtr;
//
//public:
//	template< class... Args>
//	virtual basePtr create( const Args&... args ) const = 0;
//	virtual ~Creator() = default;
//};
//
//template < class Base, class Product >
//class ProductCreator : public Creator<Base> {
//
//	typedef std::unique_ptr<Product> prodPtr;
//
//public:
//
//	template< class... Args>
//	prodPtr create( const Args&... args ) const {
//		return std::make_unique<Product>(args...);
//	}
//};
//


//template < class Base, typename T >
//class MyFactory {
//
//	typedef std::unique_ptr<MyFactory> 		facPtr;
//	typedef std::unique_ptr<Base> 			basePtr;
//	typedef std::unique_ptr<Creator<Base>> 	crePtr;
//
//private:
//
//	static std::map<T,crePtr> m_Products;
//
//public:
//
//	static facPtr get() {
//		static facPtr fac;
//		return fac;
//	}
//
//	template < T key, class Product >
//	void registerProduct( ) {
//
//		static_assert ( m_Products.find(key) == m_Products.end() ,
//				"Product is already registered in Factory" );
//
//        m_Products[key] = std::make_unique<ProductCreator<Base,Product>>();
//	}
//
//	template< class... Args >
//	basePtr create( const T & key, const Args&... args  ) const {
//
//        auto i = m_Products.find(key);
//        MYASSERT ( i != m_Products.end() , "Object is not registered!");
//
//        auto iCreator = i->second;
//        return iCreator->create(args...);
//    }
//
//
//};
//
//template < class Base >
//class MyFactoryDum {
//
//	typedef std::unique_ptr<MyFactoryDum> 		facPtr;
//	typedef std::unique_ptr<Base> 			basePtr;
//	typedef std::unique_ptr<Creator<Base>> 	crePtr;
//
//private:
//
//	static std::map<int,crePtr> m_Products;
//
//public:
//
//	static MyFactoryDum& get() {
//		static MyFactoryDum fac;
//		return fac;
//	}
//
//	template < int key, class Product >
//	void registerProduct( ) {
//
//		static_assert ( m_Products.find(key) == m_Products.end() ,
//				"Product is already registered in Factory" );
//
//        m_Products[key] = std::make_unique<ProductCreator<Base,Product>>();
//	}
//
//	template< class... Args >
//	basePtr create( const T & key, const Args&... args  ) const {
//
//        auto i = m_Products.find(key);
//        MYASSERT ( i != m_Products.end() , "Object is not registered!");
//
//        auto iCreator = i->second;
//        return iCreator->create(args...);
//    }
//
//
//};
//
//
//
////template < class Base, typename T >
////template < T key, class Product >
////void MyFactory<Base,T>::registerProduct(){
////
////	static_assert ( m_Products.find(key) == m_Products.end() ,
////					"Product is already registered in Factory" );
////
////	m_Products[key] = std::make_unique<ProductCreator<Base,Product>>();
////
////}
//
//class Dummy {
//Dummy() = default;
//virtual ~Dummy() = default;
//};
//
//class Dummy1 : public Dummy {
//
//};
//
//class Dummy2 : public Dummy {
//
//};
//
//
//typedef MyFactoryDum<Dummy> dumFac;
//
//dumFac::get();
//
//MyFactoryDum<double>::get();
//
//static MyFactoryDum<Dummy>::get();
//
//static auto fac = MyFactory<Dummy,uint8_t>::get();
//fac->registerProduct<1,Dummy1>();
//
//static MyFactory<Dummy,uint8_t>::get()->create(1);
//
//MyFactory<Dummy,uint8_t>::get()->registerProduct<1,Dummy1>();
//
//MyFactory<Dummy,uint8_t>.create(1);

//template < class Base , template<class Arg1, class Arg2> class Derived, class Arg1, class Arg2 >
//class TemplateFactory2TP {
//
//	typedef std::unique_ptr<Base>	BasePtr;
//
//private:
//
//	const std::vector<Arg1> tp1;
//	const std::vector<Arg2> tp2;
//
//public:
//
//	TemplateFactory2TP( const std::vector<Arg1> & arg1 , const std::vector<Arg2> & arg2 );
//
//	BasePtr create ( const Arg1 arg1 , const Arg2 arg2 ) {
//
//		for ( const Arg1 & p1 : tp1 ) {
//
//			for ( const Arg2 & p2 : tp2 ) {
//
//				if ( p1 == arg1 && p2 == arg2 ){
//					return BasePtr ( new Derived<arg1,arg2>() );
//				}
//			}
//
//		}
//
//		MYASSERT (false, "Dimension must be 1, 2 or 3!");
//
//
//	}
//};

#endif /* PRIMITIVES_PRIMITIVEFACTORY_H_ */

















