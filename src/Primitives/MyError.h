/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * MyError.h
 *
 *  Created on: Aug 12, 2014
 *      Author: kustepha
 */

#ifndef MYERROR_H_
#define MYERROR_H_

#define MYASSERT(condition,msg) MyError::myAssert(condition,msg,__FILE__,__LINE__)
#define MYWARN(condition,msg) MyError::myWarn(condition,msg,__FILE__,__LINE__)

#include <string>

namespace MyError
{


  __attribute__((noreturn,cold))
  void
  terminate(
      std::string msg,
      const char* file,
      const int line);

  __attribute__((cold))
  void
  warn(
      std::string msg,
      const std::string& file,
      const int line);


  __attribute__((always_inline))
  inline
  void
  myAssert(
      const bool condition,
      const std::string& msg,
      const char* file,
      const int line)
  {
    if (__builtin_expect(condition,true)) //if ( __builtin_expect_with_probability(condition,true,1.0) )
      return;
    else
      terminate(msg,file,line);
  }


  __attribute__((always_inline))
  inline
  bool
  myWarn(
      const bool condition,
      const std::string& msg,
      const std::string& file,
      const int line)
  {
    if ( __builtin_expect(condition,true) )
      return true;
    else {
        warn(msg,file,line);
        return false;
    }
  }

} /* namespace MyError */

#endif /* MYERROR_H_ */
