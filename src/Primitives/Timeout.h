/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * Timeout.h
 *
 *  Created on: Dec 8, 2017
 *      Author: kustepha
 */

#ifndef SRC_PRIMITIVES_TIMEOUT_H_
#define SRC_PRIMITIVES_TIMEOUT_H_

#include <chrono>
#include <mutex>
#include <thread>
#include <future>
#include <condition_variable>
#include <type_traits>

namespace Timeout
{
  template<typename DurationT = std::chrono::seconds>
  void throw_on_timeout(DurationT&& duration)
  {
    throw std::runtime_error(
        "timeout after waiting for "
        + std::to_string( std::chrono::duration_cast<std::chrono::duration<double>>(duration).count() )
    + "s" );
  }

  template<typename DurationT = std::chrono::seconds>
  void sleep_for(DurationT&& duration) { std::this_thread::sleep_for(duration); }

  template<class TaskT, typename DurationT = std::chrono::seconds, class TimeoutTaskT = void(*)(DurationT&&)>
  std::invoke_result_t<TaskT> wait_for(
      TaskT&& task,
      DurationT&& duration,
      TimeoutTaskT&& timeouttask = &throw_on_timeout)
      {

    std::mutex m;
    std::condition_variable cv;

    bool ready = false;

    std::thread t([&cv,&m,&duration,&ready,&timeouttask]()->void{
      std::unique_lock<std::mutex> l(m);
      if( !(cv.wait_for( l, duration, [&ready]()->bool{return ready;} )) )
        timeouttask(duration);
    });

    std::invoke_result_t<TaskT> res = task();

    // signal time-keeper thread
    {
      std::lock_guard<std::mutex> lk(m);
      ready = true;
    }
    cv.notify_one();

    t.join();

    return res;


    //    std::future<std::invoke_result_t<TaskT>> f = std::async(std::launch::async, task);
    //
    //    std::future_status status = f.wait_for(duration);
    //    if (status != std::future_status::ready)
    //      throw std::runtime_error("timeout");
    //    else
    //      return f.get();


    //    std::mutex m;
    //    std::condition_variable cv;
    //
    //    std::invoke_result_t<TaskT> res;
    //
    //    std::thread t([&cv,&m,&res,task]()->void{
    //      {
    //        std::lock_guard<std::mutex> l(m);
    //        res = task();
    //      }
    //      cv.notify_one(); });
    //
    //    t.detach();
    //
    //    {
    //      std::unique_lock<std::mutex> l(m);
    //      if( cv.wait_for( l, duration ) == std::cv_status::timeout )
    //        throw std::runtime_error("timeout");
    //      return res;
    //    }

      }

}


#endif /* SRC_PRIMITIVES_TIMEOUT_H_ */
