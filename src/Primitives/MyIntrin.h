/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * MyIntrin.h
 *
 *  Created on: Jun 20, 2014
 *      Author: Stephan Kuechlin
 */

#ifndef MYINTRIN_H_
#define MYINTRIN_H_


#include <immintrin.h>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <type_traits>

typedef double   v16df __attribute__ ((__vector_size__ (128)));
typedef double   v8df  __attribute__ ((__vector_size__ (64)));
typedef double   v4df  __attribute__ ((__vector_size__ (32)));
typedef double   v2df  __attribute__ ((__vector_size__ (16)));
typedef uint64_t v4du  __attribute__ ((__vector_size__ (32)));
typedef uint64_t v2du  __attribute__ ((__vector_size__ (16)));
typedef uint8_t  v32qu __attribute__ ((__vector_size__ (32)));
typedef int64_t  v4di  __attribute__ ((__vector_size__ (32)));
typedef float    v4sf  __attribute__ ((__vector_size__ (16)));

typedef double  d_aligned128 __attribute__ ((aligned(128)));
typedef double  d_aligned64  __attribute__ ((aligned(64)));
typedef double  d_aligned32  __attribute__ ((aligned(32)));
typedef double  d_aligned16  __attribute__ ((aligned(16)));
typedef int64_t i_aligned32  __attribute__ ((aligned(32)));

template<typename T, uint8_t len>
struct vector_type {};

template<> struct vector_type<double,  16> { using type = v16df; };
template<> struct vector_type<double,   8> { using type = v8df; };
template<> struct vector_type<double,   4> { using type = v4df; };
template<> struct vector_type<double,   2> { using type = v2df; };
template<> struct vector_type<float,    4> { using type = v4sf; };
template<> struct vector_type<uint64_t, 4> { using type = v4df; };
template<> struct vector_type<uint64_t, 2> { using type = v2du; };
template<> struct vector_type<uint8_t, 32> { using type = v32qu; };
template<> struct vector_type<int64_t,  4> { using type = v4di; };

template<typename T>
struct scalar_type {
  using type = std::remove_all_extents_t<T>;
  static_assert(sizeof(T)==sizeof(type));
};

template<> struct scalar_type<v16df> { using type = double; };
template<> struct scalar_type<v8df>  { using type = double; };
template<> struct scalar_type<v4df>  { using type = double; };
template<> struct scalar_type<v2df>  { using type = double; };
template<> struct scalar_type<v4sf>  { using type = float; };
template<> struct scalar_type<v4du>  { using type = uint64_t; };
template<> struct scalar_type<v2du>  { using type = uint64_t; };
template<> struct scalar_type<v32qu> { using type = uint8_t; };
template<> struct scalar_type<v4di>  { using type = int64_t; };
template<typename T, int len>
struct scalar_type<Eigen::Matrix<T,len,1>> { using type = T; };


template<typename T> constexpr uint8_t len() { return sizeof(T)/sizeof(typename scalar_type<T>::type); }

template<typename T>
constexpr T zero();

template<> constexpr double zero<double>() {return 0.;}
template<> constexpr v2df zero<v2df>() {return v2df{0.,0.};}
template<> constexpr v4df zero<v4df>() {return v4df{0.,0.,0.,0.};}
template<> constexpr v4du zero<v4du>() {return v4du{0,0,0,0};}

template<typename T, uint8_t n = len<T>()>
constexpr T one();

template<> constexpr double one<double,1>() {return 1.;}
template<> constexpr uint64_t one<uint64_t,1>() {return 1ull;}
template<> constexpr v4df one<v4df,1>() {return v4df{1.,0.,0.,0.};}
template<> constexpr v4df one<v4df,2>() {return v4df{1.,1.,0.,0.};}
template<> constexpr v4df one<v4df,3>() {return v4df{1.,1.,1.,0.};}
template<> constexpr v4df one<v4df,4>() {return v4df{1.,1.,1.,1.};}
template<> constexpr v4du one<v4du,1>() {return v4du{1,0,0,0};}
template<> constexpr v4du one<v4du,2>() {return v4du{1,1,0,0};}
template<> constexpr v4du one<v4du,3>() {return v4du{1,1,1,0};}
template<> constexpr v4du one<v4du,4>() {return v4du{1,1,1,1};}

//template<typename T>
//constexpr auto& to_eigen(T& v) { return *(Eigen::Matrix<typename scalar_type<T>::type,len<T>(),1>*)(&v); }
//
//template<typename T>
//constexpr const auto& to_eigen(const T& v) { return *(const Eigen::Matrix<typename scalar_type<T>::type,len<T>(),1>*)(&v); }
//
//template<typename T, int len>
//constexpr auto& from_eigen(Eigen::Matrix<T,len,1>& v) { return *(typename vector_type<T,uint8_t(len)>::type*)(&v); }
//
//template<typename T, int len>
//constexpr const auto& from_eigen(const Eigen::Matrix<T,len,1>& v) { return *(const typename vector_type<T,uint8_t(len)>::type*)(&v); }

template<typename T>
__attribute__((__always_inline__))
inline constexpr auto to_eigen(const T& v) {
  using ST = typename scalar_type<T>::type;
  typename Eigen::Matrix<ST,len<T>(),1> r;
  __builtin_memcpy(r.data(),&v,len<T>()*sizeof(ST));
  return r;
}

template<typename T, int len>
__attribute__((__always_inline__))
inline constexpr auto from_eigen(const Eigen::Matrix<T,len,1>& v) {
  typename vector_type<T,uint8_t(len)>::type r;
  __builtin_memcpy(&r,&v,len*sizeof(T));
  return r;
}


template<uint8_t n>
constexpr auto bits(uint8_t v) {
  typename vector_type<int64_t,n>::type r;
  for (int8_t i = 0; i < n; ++i, v>>=1) r[i] = v&1;
  return r;
}

template<uint8_t n>
constexpr v4du active4du() {
  switch (n%4) {
    case 0: return v4du{0,0,0,0}; break;
    case 1: return v4du{1,0,0,0}; break;
    case 2: return v4du{1,1,0,0}; break;
    case 3: return v4du{1,1,1,0}; break;
  }
  __builtin_unreachable();
  return zero<v4du>();
}

template<uint8_t n>
constexpr v4df active4df() {
  switch (n%4) {
    case 0: return v4df{0.,0.,0.,0.}; break;
    case 1: return v4df{1.,0.,0.,0.}; break;
    case 2: return v4df{1.,1.,0.,0.}; break;
    case 3: return v4df{1.,1.,1.,0.}; break;
  }
  __builtin_unreachable();
  return zero<v4df>();
}

template<uint8_t n>
constexpr v4du inactive4du() {
  return v4du{1,1,1,1} - active4du<n>();  }

template<uint8_t n>
constexpr v4df inactive4df() {
  return v4df{1.,1.,1.,1.} - active4df<n>();  }

inline std::ostream&
operator<< (std::ostream& s, const v4df& v){
  return s << v[0] << " " << v[1] << " " << v[2] << " " << v[3]; }

//inline std::ostream&
//operator<< (std::ostream& s, const v2df& v){
//  return s << v[0] << " " << v[1]; }

inline std::ostream&
operator<< (std::ostream& s, const v4di& v){
  return s << v[0] << " " << v[1] << " " << v[2] << " " << v[3];}

inline std::ostream&
operator<< (std::ostream& s, const v4du& v){
  return s << v[0] << " " << v[1] << " " << v[2] << " " << v[3];}

inline
std::string to_string(const v4df& v) {
  return std::to_string((double)(v[0])) + " "
      + std::to_string((double)(v[1])) + " "
      + std::to_string((double)(v[2])) + " "
      + std::to_string((double)(v[3])); }

template<size_t sz>
struct reg_t {};

template<>
struct reg_t<16> {
  static constexpr size_t size = 16;
  typedef __m128i type;
  type reg;
  void load(void* addr) { reg = _mm_stream_load_si128((type*)addr); }
  void store(void* addr) const { _mm_stream_si128((type*)addr,reg); }
};

#if defined(__AVX__)
template<>
struct reg_t<32> {
  static constexpr size_t size = 32;
  typedef __m256i type;
  type reg;
  void load(void* addr) {
#if defined(__AVX2__)
    reg = _mm256_stream_load_si256((type*)addr);
#else
    reg = _mm256_load_si256((type*)addr);
#endif
  }
  void store(void* addr) const { _mm256_stream_si256((type*)addr,reg); }
};
#endif

#if defined(__AVX512F__)
template<>
struct reg_t<64> {
  static constexpr size_t size = 64;
  typedef __m512i type;
  type reg;
  void load(void* addr) { reg =
      //       _mm512_load_si512((type*)addr);
      _mm512_stream_load_si512((type*)addr);
  }
  void store(void* addr) const {
    //_mm512_store_si512((type*)addr,reg);
    _mm512_stream_si512((type*)addr,reg);
  }
};

#endif


namespace MyIntrin {


  __attribute__((__always_inline__))
          inline v4du gather(const v4du idx, const uint64_t* const src)
  {
    __m256i res = _mm256_setzero_si256 ();

#if defined( __AVX__ )

    size_t r;

    r = idx[0];
    res = _mm256_insert_epi64(res,src[r],0);
    r = idx[1];
    res = _mm256_insert_epi64(res,src[r],1);
    r = idx[2];
    res = _mm256_insert_epi64(res,src[r],2);
    r = idx[3];
    res = _mm256_insert_epi64(res,src[r],3);

#else
    res[ 0 ] = src[ idx[0] ];
    res[ 1 ] = src[ idx[1] ];
    res[ 2 ] = src[ idx[2] ];
    res[ 3 ] = src[ idx[3] ];
#endif

    return (v4du)res;
  }


  __attribute__((__always_inline__))
  inline void scatter(const v4du& __restrict__ idx, const v4du& __restrict__ src, uint64_t* const dst )
  {

#if defined( __AVX__ )

    __m256i tmp0, tmp1, tmp2, tmp3;
    size_t r;

    __asm__
    __volatile__ // tell optimizer not to ignore this
    (// GAS syntax: mnemonic source, destination

        "vmovq %x[idx], %[r] \n"                 // r <- idx[0]
        "vmovq %x[src], (%[dst],%[r],8) \n"      // dst[r] <- src[0]
        "vpalignr $8,%x[idx],%x[idx],%x[tmp0] \n"   // tmp0 <-idx[0],idx[1]
        "vpalignr $8,%x[src],%x[src],%x[tmp1] \n"   // tmp1 <-src[0],src[1]
        "vmovq %x[tmp0], %[r] \n"                // r <- idx[1]
        "vmovq %x[tmp1], (%[dst],%[r],8) \n"     // dst[r] <- src[1]

        "vextractf128 $1,%t[idx],%x[tmp0] \n"       // tmp0 <- idx[3],idx[2]
        "vextractf128 $1,%t[src],%x[tmp1] \n"       // tmp1 <- src[3],src[2]

        "vmovq %x[tmp0], %[r] \n"                // r <- idx[2]
        "vmovq %x[tmp1], (%[dst],%[r],8) \n"     // dst[r] <- src[2]
        "vpalignr $8,%x[tmp0],%x[tmp0],%x[tmp2] \n" // tmp2 <-idx[2],idx[3]
        "vpalignr $8,%x[tmp1],%x[tmp1],%x[tmp3] \n" // tmp3 <-src[2],src[3]
        "vmovq %x[tmp2], %[r] \n"                // r <- idx[3]
        "vmovq %x[tmp3], (%[dst],%[r],8) \n"     // dst[r] <- src[3]
        : [r]"=&r"(r),
          [tmp0]"=&x"(tmp0),[tmp1]"=&x"(tmp1),
          [tmp2]"=&x"(tmp2),[tmp3]"=&x"(tmp3) // temporaries
          : [idx]"x"(idx),[src]"x"(src), [dst]"r"(dst)      //input
            : "memory" // clobbers
    );

#else
    dst[ idx[0] ] = src[0];
    dst[ idx[1] ] = src[1];
    dst[ idx[2] ] = src[2];
    dst[ idx[3] ] = src[3];
#endif
  }


  __attribute__((__always_inline__))
  inline v4df cvt_du2df(const v4du& v)
  {
#if defined(__AVX512VL__) && defined(__AVX512DQ__)
    return _mm256_cvtepu64_pd ( (__m256i)v );
#else
    return v4df{
      static_cast<double>(v[0]),
          static_cast<double>(v[1]),
          static_cast<double>(v[2]),
          static_cast<double>(v[3])
    };
#endif
  }

  __attribute__((__always_inline__))
  inline v4du cvt_df2du_truncate(const v4df& v)
  {
#if defined(__AVX512VL__) && defined(__AVX512DQ__)
    return (v4du)_mm256_cvtpd_epu64 (
        _mm256_round_pd( (__m256d)v, (_MM_FROUND_TO_ZERO |_MM_FROUND_NO_EXC) ) );
#else
    return v4du{
      static_cast<uint64_t>(v[0]),
          static_cast<uint64_t>(v[1]),
          static_cast<uint64_t>(v[2]),
          static_cast<uint64_t>(v[3])
    };
#endif
  }

  __attribute__((__always_inline__))
  inline v4df set(const double v0, const double v1, const double v2, const double v3)
  {
#ifdef __AVX__
    return _mm256_set_pd(v3,v2,v1,v0);
#else
    return v4df{v0,v1,v2,v3};
#endif
  }

  template<typename VT, typename FT, typename... AT>
  __attribute__((__always_inline__))
  inline VT forall(const FT&& f, const VT v, AT... args) {
    VT ret;
    for (int i = 0; i < len<VT>(); ++i)
      ret[i] = f((double)v[i],args...);
    return ret;
  }

  template<typename VT, typename FT, typename... AT>
  __attribute__((__always_inline__))
  void transform(const FT&& f, const VT v, AT... args) {
    for (int i = 0; i < len<VT>(); ++i)
      f((double)v[i],args...);
  }

  template<typename VT, uint8_t n = len<VT>()>
  __attribute__((__always_inline__))
  inline bool all(const VT& v) {
    static_assert(std::is_same<typename scalar_type<VT>::type,int64_t>::value);
    static_assert(n<=len<VT>());
    int64_t ret = -1;
    for (uint8_t i =0; i < n; ++i)
      ret &= v[i];
    return ret!=0;
  }

  __attribute__((__always_inline__))
  inline v4df sqrt4(const v4df v) {
    return _mm256_sqrt_pd(v); }


  __attribute__((__always_inline__))
  inline v4df pow4(const v4df base, const double power)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = pow((double)base[i],power);
    return res;
  }

  __attribute__((__always_inline__))
  inline v4df exp4(const v4df v)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = exp((double)v[i]);
    return res;
  }

  __attribute__((__always_inline__))
  inline v4df expm14(const v4df v)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = expm1((double)v[i]);
    return res;
  }

  __attribute__((__always_inline__))
  inline v4df log4(const v4df v)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = log((double)v[i]);
    return res;
  }


  __attribute__((__always_inline__))
  inline v4df cos4(const v4df v)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = cos((double)v[i]);
    return res;
  }

  __attribute__((__always_inline__))
  inline v4df sin4(const v4df v)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = sin((double)v[i]);
    return res;
  }

  __attribute__((__always_inline__))
  inline v4df tanh4(const v4df v)
  {
    v4df res;
    for (int i = 0; i < 4; ++i)
      res[i] = tanh((double)v[i]);
    return res;
  }


  template<typename VT>
  __attribute__((__always_inline__))
  inline
  VT round(
      VT v,
      const double tol,
      const double invtol)
  {
    v *= invtol;

    v = forall([](double x){return ::round(x);},v);

    return v*tol;
  }

  __attribute__((__always_inline__))
  inline v4df set(const double v)
  {
#ifdef __AVX__
    return _mm256_set1_pd(v);
#else
    return v4df{v,v,v,v};
#endif
  }

  __attribute__((__always_inline__))
  inline v4df load4du(const double *const P)
  {
#ifdef __AVX__
    return _mm256_loadu_pd(P);
#else
    return v4df{P[0], P[1], P[2], P[3]};
#endif
  }

  __attribute__((__always_inline__))
  inline v4df load4d(const double *const P)
  {
#ifdef __AVX__
    return _mm256_load_pd(P);
#else
    return v4df{P[0], P[1], P[2], P[3]};
#endif
  }

#ifdef EIGEN_WORLD_VERSION
  __attribute__((__always_inline__))
  inline v4df load4d(const Eigen::Vector4d& v)
  {
#ifdef __AVX__
    return _mm256_load_pd(v.data());
#else
    return *(v4df*)v.data();
#endif
  }
#endif

  __attribute__((__always_inline__))
  inline v4du load4u(const uint64_t *const P)
  {
#ifdef __AVX__
    return (v4du)_mm256_load_si256((__m256i*)P);
#else
    return *(v4du*)P;
#endif
  }

  __attribute__((__always_inline__))
  inline v4di load4i(const int64_t *const P)
  {
#ifdef __AVX__
    return (v4di)_mm256_load_si256((__m256i*)P);
#else
    return *(v4di*)P;
#endif
  }

  __attribute__((__always_inline__))
  inline void store4d(double *const P, const v4df& v)
  {
#ifdef __AVX__
    _mm256_store_pd(P,v);
#else
    *P = v;
#endif
  }

  __attribute__((__always_inline__))
  inline void stream4d(double *const P, const v4df& v)
  {
#ifdef __AVX__
    _mm256_stream_pd(P,v);
#else
    *P = v;
#endif
  }

  __attribute__((__always_inline__))
  inline void stream4i(int64_t *const P, const v4di& v)
  {
#ifdef __AVX__
    _mm256_stream_si256((__m256i*)P,(__m256i)v);
#else
    *P = v;
#endif
  }

  __attribute__((__always_inline__))
  inline v4di streamload4i(const int64_t *const P)
  {
#ifdef __AVX2__
    return (v4di)_mm256_stream_load_si256((__m256i*)P);
#else
    return load4i(P);
#endif
  }

  template<char sz>
  __attribute__((__always_inline__))
  inline
  void copy(void* dst, const void* src)
  {
    typedef char vec __attribute__ ((__vector_size__ (sz)));
    *(vec*)dst = *(vec*)src;
  }


  template<uint8_t i>
  __attribute__((__always_inline__))
  inline constexpr double get(const v4df& v) {
    return v[i]; }

  __attribute__((__always_inline__))
  inline constexpr double get(const v4df& v, const uint i) {
    return v[i]; }


  inline bool isfinite(const v4df& v) {
    return  std::isfinite(get(v,0)) && std::isfinite(get(v,1)) &&
        std::isfinite(get(v,2)) && std::isfinite(get(v,3)); }



  __attribute__((__always_inline__))
  inline v4df set(const double* const P)
  {
#ifdef __AVX__
    return  _mm256_broadcast_sd(P);
#else
    return v4df{*P, *P, *P, *P};
#endif
  }

  __attribute__((__always_inline__))
  inline v4df set3(const double P)
  {
    return set(P,P,P,0.0);
  }

  __attribute__((__always_inline__))
  inline v4df set3(const double* const P)
  {
    return set3(*P);
  }

  /*
  __attribute__((__always_inline__))
  inline constexpr double shsum( const v2df& v ) { return v[0] + v[1]; }
  */

  __attribute__((__always_inline__))
  inline constexpr double shsum( const v4df& v ) { return v[0] + v[1] + v[2] + v[3]; }

  __attribute__((__always_inline__))
  inline constexpr double shsum3( const v4df& v ) { return v[0] + v[1] + v[2]; }

  /*
  inline constexpr v2df hsum( const v2df& v ) {
    return v + __builtin_shuffle(v,v2du{1,0}); }
    */

  inline constexpr v4df hsum( const v4df& v ) {
    v4df s0 = v + __builtin_shuffle(v,v4du{1,0,3,2});
    return s0 + __builtin_shuffle(s0,v4du{2,3,0,1}); }

  inline constexpr v4df hsum2( const v4df& v ) {
    v4df v3 = __builtin_shuffle(v,zero<v4df>(),v4du{0,1,4,5});
    return v3 + __builtin_shuffle(v3,zero<v4df>(),v4du{1,0,2,3}); }

  inline constexpr v4df hsum3( const v4df& v ) {
    return hsum(v*active4df<3>())*active4df<3>(); }

  inline constexpr void bcast3(
      const v4df& src,
      v4df& d0,
      v4df& d1,
      v4df& d2   )
  {
    constexpr v4du mask0 = {0,0,0,0};
    constexpr v4du mask1 = {1,1,1,1};
    constexpr v4du mask2 = {2,2,2,2};

    d0 = __builtin_shuffle(src,mask0);
    d1 = __builtin_shuffle(src,mask1);
    d2 = __builtin_shuffle(src,mask2);
  }

  inline constexpr void bcast4(
      const v4df& src,
      v4df& d0,
      v4df& d1,
      v4df& d2,
      v4df& d3 )
  {
    constexpr v4du mask0 = {0,0,0,0};
    constexpr v4du mask1 = {1,1,1,1};
    constexpr v4du mask2 = {2,2,2,2};
    constexpr v4du mask3 = {3,3,3,3};

    d0 = __builtin_shuffle(src,mask0);
    d1 = __builtin_shuffle(src,mask1);
    d2 = __builtin_shuffle(src,mask2);
    d3 = __builtin_shuffle(src,mask3);
  }


  inline constexpr v4df symmm3x3_times_v3(
      const v4df& m0j,
      const v4df& m1j,
      const v4df& m2j,
      const v4df& v)
  {
    v4df v0(zero<v4df>());
    v4df v1(zero<v4df>());
    v4df v2(zero<v4df>());
    bcast3(v,v0,v1,v2);
    return m0j*v0 + m1j*v1 + m2j*v2;
  }

  inline constexpr v4df symmm3x3_times_v3(
      const double m[6],
      const v4df& v)
  {
    v4df v0(zero<v4df>());
    v4df v1(zero<v4df>());
    v4df v2(zero<v4df>());
    bcast3(v,v0,v1,v2);
    return \
        v4df{m[0],m[1],m[2],0.}*v0 +
        v4df{m[1],m[3],m[4],0.}*v1 +
        v4df{m[2],m[4],m[5],0.}*v2;
  }

  inline constexpr v4df symmm3x3_mkl_vkvl(
      const double m[6],
      const v4df& v)
  {
    return hsum3(
        v*v*v4df{m[0],m[3],m[5],0.} +
        2.*v4df{m[1],m[2],m[4],0.} *
        v4df{v[1],v[0],v[2],v[3]} *
        v4df{v[0],v[2],v[1],v[3]});
  }


  inline constexpr v4df symmm4x4_times_v4(
      const v4df& m0j,
      const v4df& m1j,
      const v4df& m2j,
      const v4df& m3j,
      const v4df& v)
  {
    v4df v0(zero<v4df>());
    v4df v1(zero<v4df>());
    v4df v2(zero<v4df>());
    v4df v3(zero<v4df>());
    bcast4(v,v0,v1,v2,v3);
    return m0j*v0 + m1j*v1 + m2j*v2 + m3j*v3;
  }


  /*
  __attribute__((__always_inline__))
  inline constexpr double sdot(
      const v2df& va,
      const v2df& vb )
  {
    return shsum(va*vb);
  }
  */

  __attribute__((__always_inline__))
  inline constexpr double sdot(
      const v4df& va,
      const v4df& vb )
  {
    return shsum(va*vb);
  }

  __attribute__((__always_inline__))
  inline constexpr double sdot3(
      const v4df& va,
      const v4df& vb )
  {
    return shsum3(va*vb);
  }

  /*
  __attribute__((__always_inline__))
  inline constexpr v2df dot(const v2df& va, const v2df& vb )
  {
    return hsum( va * vb );
  }
  */

  __attribute__((__always_inline__))
  inline constexpr v4df dot(const v4df& va, const v4df& vb )
  {
    return hsum( va * vb );
  }

  __attribute__((__always_inline__))
  inline constexpr v4df dot2(const v4df& va, const v4df& vb )
  {
    return hsum2( va * vb );
  }


  __attribute__((__always_inline__))
  inline constexpr v4df dot3(const v4df& va, const v4df& vb )
  {
    return hsum3( va * vb );
  }

  inline v4df cross3(const v4df& va, const v4df& vb)
  {
    static constexpr v4du mask0 = {1,2,0,3};
    static constexpr v4du mask1 = {2,0,1,3};
    v4df a0 = __builtin_shuffle(va,mask0);
    v4df a1 = __builtin_shuffle(va,mask1);
    v4df b0 = __builtin_shuffle(vb,mask0);
    v4df b1 = __builtin_shuffle(vb,mask1);
    return a0*b0 - a1*b1;
  }

  /*
  __attribute__((__always_inline__))
  inline constexpr v2df abssq(const v2df& v) { return dot(v,v); }
*/

  __attribute__((__always_inline__))
  inline constexpr v4df abssq(const v4df& v) { return dot(v,v); }

  __attribute__((__always_inline__))
  inline constexpr v4df abssq3(const v4df& v) { return dot3(v,v); }

  /*
  __attribute__((__always_inline__))
  inline constexpr double sabssq(const v2df& v) { return sdot(v,v); }
*/

  template<typename T,uint N = sizeof(T)/sizeof(double)>
  __attribute__((__always_inline__))
  inline constexpr double sabssq(const T& v) {
    double ans = v[0]*v[0];
    for (uint i = 1; i < N; ++i) ans += v[i]*v[i];
    return ans;
  }

  __attribute__((__always_inline__))
  inline constexpr double sabssq(const v4df& v) { return sdot(v,v); }

  __attribute__((__always_inline__))
  inline constexpr double sabssq3(const v4df& v) { return sdot3(v,v); }

  __attribute__((__always_inline__))
  inline double snorm(const v4df& v) { return sqrt(sabssq(v)); }

  inline v4df prod(v4df v)
  {
    v *= __builtin_shuffle(v,v4du{1,0,3,2});
    return v *= __builtin_shuffle(v,v4du{2,3,0,1});
  }

  template <uint8_t n>
  static constexpr double sprod(const v4df& v) {
    double r = 1.;
    for (uint8_t i = 0; i < n; ++i)
      r *= v[i];
    return r;
  }

  template <uint8_t n>
  static constexpr uint64_t sprod(const v4du& v) {
    uint64_t r = 1.;
    for (uint8_t i = 0; i < n; ++i)
      r *= v[i];
    return r;
  }

  __attribute__((__always_inline__))
  inline v4df rot(
      const v4df& v )
  {
    static constexpr v4du mask = {3,0,1,2};
    return __builtin_shuffle(v,mask);
  }


  inline bool is_approx(const v4df& a, const v4df& b, const double reltol) {
    constexpr struct { double constexpr operator() (const double v) const { return v != 0.0 ? v : 1.0; } } avoid0;
    return reltol > (2.0*sabssq(a - b) / avoid0( sabssq(a + b) ));
  }

  __attribute__((__always_inline__))
  inline v4di max(const v4di& a, const v4di& b)
  {
    return a ^ ( (a^b) & (b>a) );
  }

  __attribute__((__always_inline__))
  inline v4df max(const v4df& a, const v4df& b)
  {
#ifdef __AVX__
    return (v4df)_mm256_max_pd((__m256d)a,(__m256d)b);
#else
    return a > b ? a : b;
#endif
  }

  __attribute__((__always_inline__))
  inline v4di hmax(const v4di& a)
  {
    //#ifdef __AVX__
    //  v4di b = _mm256_permute2f128_si256(a,a,1); // b = a[2] a[3] a[0] a[1]
    //  v4di m1 = max(a,b); // m1 = max(a[0],a[2]) max(a[1],a[3]) max(a[2],a[0]) max(a[1],a[3])
    //  v4di m2 = (v4di)_mm256_permute_pd((v4df)m1,5); // m2 = max(a[1],a[3]), max(a[0],a[2]), max(a[1],a[3]), max(a[2],a[0])
    //  return max(m2,m1);
    //#else
    v4di b = __builtin_shuffle(a,v4di{2,3,0,1});
    v4di m1 = max(a,b);
    v4di m2 = __builtin_shuffle(m1,v4di{1,0,3,2});
    return max(m2,m1);
    //#endif
  }

  template<typename T, uint8_t n = sizeof(T)/sizeof(typename scalar_type<T>::type)>
  __attribute__((__always_inline__))
  inline auto shmin(const T& v) {
    auto m = v[0];
    for (uint8_t i = 1; i < n; ++i)
      m = std::min(m,v[i]);
    return m;
  }

  template<typename T>
  __attribute__((__always_inline__))
  inline auto shmax(const T& v) {
    auto m = v[0];
    for (uint8_t i = 1; i < sizeof(T)/sizeof(m); ++i)
      m = std::max(m,v[i]);
    return m;
  }


#ifdef __AVX2__
  inline
  __m256i byte_rank4(const __m256i v)
  {
    __m256i res0 = __m256i{-3,-2,-1,0};
    __m256i cmp0 = v == __builtin_shuffle(v,__m256i{0,0,1,2});
    __m256i cmp1 = v == __builtin_shuffle(v,__m256i{0,1,0,1});
    __m256i cmp2 = v == __builtin_shuffle(v,__m256i{0,1,2,0});
    __m256i res1 = _mm256_sub_epi64(res0,cmp0);
    __m256i res2 = _mm256_sub_epi64(res1,cmp1);
    __m256i res3 = _mm256_sub_epi64(res2,cmp2);
    return res3;
  }

  inline
  v32qu byte_rank32(const v32qu v)
  {
    const __m256i corr = _mm256_set_epi8(
        -1, -2, -3, -4, -5, -6, -7,
        -8, -9,-10,-11,-12,-13,-14,-15,
        -16,-17,-18,-19,-20,-21,-22,-23,
        -24,-25,-26,-27,-28,-29,-30,-31,-32);

    const __m256i zero = _mm256_setzero_si256();

    __m256i vv,res1,res2,vrot1,vrot2,cmp1,cmp2;

    vv = (__m256i)v;

    cmp1 = _mm256_cmpeq_epi8((__m256i)v,zero);
    res1 = _mm256_blendv_epi8(zero,corr,cmp1);
    res2 = zero;

    for (int8_t i = 0; i < 16; i+=2)
      {
        vrot1 =  _mm256_alignr_epi8(vv,zero,i);
        vrot2 =  _mm256_alignr_epi8(vv,zero,i+1);
        cmp1 = _mm256_cmpeq_epi8((__m256i)v,vrot1);
        cmp2 = _mm256_cmpeq_epi8((__m256i)v,vrot2);
        res1 = _mm256_sub_epi8(res1,cmp1);
        res2 = _mm256_sub_epi8(res2,cmp2);
      }

    vv = _mm256_permute2x128_si256(vv,zero,0x0F);

    for (int8_t i = 0; i < 16; i+=2)
      {
        vrot1 =  _mm256_alignr_epi8(vv,vv,i);
        vrot2 =  _mm256_alignr_epi8(vv,vv,i+1);
        cmp1 = _mm256_cmpeq_epi8((__m256i)v,vrot1);
        cmp2 = _mm256_cmpeq_epi8((__m256i)v,vrot2);
        res1 = _mm256_sub_epi8(res1,cmp1);
        res2 = _mm256_sub_epi8(res2,cmp2);
      }

    return (v32qu)_mm256_add_epi8(res1,res2);
  }

  inline
  v32qu byte_occurance32(const v32qu v)
  {

    __m256i vv,res1,res2,vrot1,vrot2,cmp1,cmp2;

    vv = (__m256i)v;

    res1 = _mm256_setzero_si256();
    res2 = _mm256_setzero_si256();

    for (int8_t i = 0; i < 16; i+=2)
      {
        vrot1 =  _mm256_alignr_epi8(vv,vv,i);
        vrot2 =  _mm256_alignr_epi8(vv,vv,i+1);
        cmp1 = _mm256_cmpeq_epi8((__m256i)v,vrot1);
        cmp2 = _mm256_cmpeq_epi8((__m256i)v,vrot2);
        res1 = _mm256_sub_epi8(res1,cmp1);
        res2 = _mm256_sub_epi8(res2,cmp2);
      }

    vv = _mm256_permute2x128_si256(vv,vv,1);

    for (int8_t i = 0; i < 16; i+=2)
      {
        vrot1 =  _mm256_alignr_epi8(vv,vv,i);
        vrot2 =  _mm256_alignr_epi8(vv,vv,i+1);
        cmp1 = _mm256_cmpeq_epi8((__m256i)v,vrot1);
        cmp2 = _mm256_cmpeq_epi8((__m256i)v,vrot2);
        res1 = _mm256_sub_epi8(res1,cmp1);
        res2 = _mm256_sub_epi8(res2,cmp2);
      }

    return (v32qu)_mm256_add_epi8(res1,res2);
  }


  inline
  __m256i gather4x4bit(__m256i v00_03, __m256i v04_07, __m256i v08_11, __m256i v12_15,
                       __m256i x)
  {

    __m256i zero = _mm256_setzero_si256();
    __m256i cmpgt03 = x > 3;
    __m256i cmpgt07 = x > 7;
    __m256i cmpgt11 = x > 11;

    __m256i res0 = __builtin_shuffle(v00_03,x);
    res0 = _mm256_blendv_epi8(res0,zero,cmpgt03);
    __m256i res1 = __builtin_shuffle(v04_07,x);
    res1 = _mm256_blendv_epi8(zero,res1,cmpgt03 & ~cmpgt07);
    __m256i res2 = __builtin_shuffle(v08_11,x);
    res2 = _mm256_blendv_epi8(zero,res2,cmpgt07 & ~cmpgt11);
    __m256i res3 = __builtin_shuffle(v12_15,x);
    res3 = _mm256_blendv_epi8(zero,res3,cmpgt11);

    return res0 ^ res1 ^ res2 ^ res3;
  }

  inline
  void scatter_4xpacked16x4bit(const __m256i x,
                               __m256i& v00_03, __m256i& v04_07, __m256i& v08_11, __m256i& v12_15)
  {
    const __m256i mask = __m256i{0xF,0xF0,0xF00,0xF000};
    const __m256i shft = __m256i{0,4,8,12};

    __m256i x1 = _mm256_alignr_epi8 (x, x, 8);
    __m256i x2 = _mm256_add_epi64(x,x1);
    __m256i x3 = _mm256_permute2x128_si256(x2,x2,1);
    __m256i h_packed16x4bit = _mm256_add_epi64(x2,x3); // hadd
    // each lane now contains the same information: 16x4bit counts

    __m256i h_packed12x4bit = _mm256_srli_epi64(h_packed16x4bit,16);
    __m256i h_packed8x4bit  = _mm256_srli_epi64(h_packed16x4bit,32);
    __m256i h_packed4x4bit  = _mm256_srli_epi64(h_packed16x4bit,48);

    v00_03  = _mm256_and_si256(h_packed16x4bit,mask);
    v04_07  = _mm256_and_si256(h_packed12x4bit,mask);
    v08_11  = _mm256_and_si256(h_packed8x4bit, mask);
    v12_15  = _mm256_and_si256(h_packed4x4bit, mask);

    v00_03  = _mm256_srlv_epi64 (v00_03, shft);
    v04_07  = _mm256_srlv_epi64 (v04_07, shft);
    v08_11  = _mm256_srlv_epi64 (v08_11, shft);
    v12_15  = _mm256_srlv_epi64 (v12_15, shft);

  }

  inline
  void hist4bit(const uint64_t* const input, const uint64_t length, const int8_t shift, uint64_t* const hist) {

    if (!length) return;

    constexpr int width = 4;
    constexpr int interleave = 3;
    constexpr int stride = width*interleave;

    __m256i h0,h4,h8,h12;

    const __m256i one  = _mm256_set1_epi64x(1);
    const __m256i mask = _mm256_set1_epi64x(0xF);
    const __m128i scount = _mm_set_epi64x(0,shift);

    const uint64_t misalign = (uint64_t)input & 31;
    const uint64_t prolog = std::min(
        length,
        misalign ? (4 - (misalign >> 3)) : 0 );
    const uint64_t epilog = prolog + stride*((length-prolog)/stride);

    // prolog
    uint64_t i = 0;
    for (; i < prolog; ++i)
      ++hist[ (input[i] >> shift) & 0xF ];

    // main SIMD loop

    __m256i hist0_3   = _mm256_load_si256((const __m256i*)&hist[0]);
    __m256i hist4_7   = _mm256_load_si256((const __m256i*)&hist[4]);
    __m256i hist8_11  = _mm256_load_si256((const __m256i*)&hist[8]);
    __m256i hist12_15 = _mm256_load_si256((const __m256i*)&hist[12]);

    for (; i < epilog; i+=stride) {

        __m256i in0 = _mm256_load_si256((const __m256i*)&input[i]);
        __m256i in1 = _mm256_load_si256((const __m256i*)&input[i+width]);
        __m256i in2 = _mm256_load_si256((const __m256i*)&input[i+2*width]);

        in0 = _mm256_srl_epi64(in0,scount);
        in1 = _mm256_srl_epi64(in1,scount);
        in2 = _mm256_srl_epi64(in2,scount);

        in0 = _mm256_and_si256(in0,mask);
        in1 = _mm256_and_si256(in1,mask);
        in2 = _mm256_and_si256(in2,mask);

        in0 = _mm256_slli_epi64(in0,2); // "x4"
        in1 = _mm256_slli_epi64(in1,2); // "x4"
        in2 = _mm256_slli_epi64(in2,2); // "x4"

        __m256i v0 = _mm256_sllv_epi64(one,in0); // "scatter"
        __m256i v1 = _mm256_sllv_epi64(one,in1); // "scatter"
        __m256i v2 = _mm256_sllv_epi64(one,in2); // "scatter"

        __m256i v3 = _mm256_add_epi64(v0,v1);
        __m256i v4 = _mm256_add_epi64(v2,v3);

        scatter_4xpacked16x4bit(v4, h0,h4,h8,h12);

        hist0_3   = _mm256_add_epi64(hist0_3,h0);
        hist4_7   = _mm256_add_epi64(hist4_7,h4);
        hist8_11  = _mm256_add_epi64(hist8_11,h8);
        hist12_15 = _mm256_add_epi64(hist12_15,h12);

    }

    _mm256_store_si256((__m256i*)&hist[0],hist0_3);
    _mm256_store_si256((__m256i*)&hist[4],hist4_7);
    _mm256_store_si256((__m256i*)&hist[8],hist8_11);
    _mm256_store_si256((__m256i*)&hist[12],hist12_15);

    // epiloge
    for (; i < length; ++i)
      ++hist[ (input[i] >> shift) & 0xF ];

  }

  template<typename SFunT, typename VFunT>
  inline
  void histfun4bit(
      const uint64_t* const input,
      const uint64_t length,
      const int8_t shift,
      uint64_t* const hist,
      const SFunT sfun,
      const VFunT vfun) {

    if (!length) return;

    constexpr int width = 4;
    constexpr int interleave = 1;
    constexpr int stride = width*interleave;

    __m256i h0,h4,h8,h12;

    const __m256i one  = _mm256_set1_epi64x(1);
    const __m256i mask = _mm256_set1_epi64x(0xF);
    const __m128i scount = _mm_set_epi64x(0,shift);

    const uint64_t misalign = (uint64_t)input & 31;
    const uint64_t prolog = std::min(
        length,
        misalign ? (4 - (misalign >> 3)) : 0 );
    const uint64_t epilog = prolog + stride*((length-prolog)/stride);

    // prolog
    uint64_t i = 0;
    for (; i < prolog; ++i)
      {
        const uint8_t v = (input[i] >> shift) & 0xF;
        sfun(i,hist[v]);
        ++hist[v];
      }

    // main SIMD loop

    __m256i hist0_3   = _mm256_load_si256((const __m256i*)&hist[0]);
    __m256i hist4_7   = _mm256_load_si256((const __m256i*)&hist[4]);
    __m256i hist8_11  = _mm256_load_si256((const __m256i*)&hist[8]);
    __m256i hist12_15 = _mm256_load_si256((const __m256i*)&hist[12]);

    for (; i < epilog; i+=stride) {

        __m256i in0 = _mm256_load_si256((const __m256i*)&input[i]);

        in0 = _mm256_srl_epi64(in0,scount);
        in0 = _mm256_and_si256(in0,mask);

        __m256i h = gather4x4bit(hist0_3,hist4_7,hist8_11,hist12_15,in0);
        h = _mm256_add_epi64(h,byte_rank4(in0));
        vfun(i,h);

        in0 = _mm256_slli_epi64(in0,2); // "x4"

        __m256i v0 = _mm256_sllv_epi64(one,in0); // "scatter"

        scatter_4xpacked16x4bit(v0, h0,h4,h8,h12);

        hist0_3   = _mm256_add_epi64(hist0_3,h0);
        hist4_7   = _mm256_add_epi64(hist4_7,h4);
        hist8_11  = _mm256_add_epi64(hist8_11,h8);
        hist12_15 = _mm256_add_epi64(hist12_15,h12);

    }

    _mm256_store_si256((__m256i*)&hist[0],hist0_3);
    _mm256_store_si256((__m256i*)&hist[4],hist4_7);
    _mm256_store_si256((__m256i*)&hist[8],hist8_11);
    _mm256_store_si256((__m256i*)&hist[12],hist12_15);

    // epiloge
    for (; i < length; ++i)
      {
        const uint8_t v = (input[i] >> shift) & 0xF;
        sfun(i,hist[v]);
        ++hist[v];
      }

  }

#endif


}// namespace MyIntrin



#endif /* MYINTRIN_H_ */
