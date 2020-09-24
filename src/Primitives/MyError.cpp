/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * MyError.cpp
 *
 *  Created on: Aug 12, 2014
 *      Author: kustepha
 */

#include "MyError.h"

#include <iostream>
#include <ctime>
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>
#include <stacktrace.h>

#ifdef DUMP_CORE

#include <google/coredumper.h>

#endif

#include "Timeout.h"

namespace MyError
{


#ifdef DUMP_CURE
  std::string dump() {

    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    char* scratchDir = getenv("SCRATCH");

    time_t t = std::time(0);
    struct tm * now = localtime( & t );

    std::string corefilename = "";
    if (scratchDir != nullptr)
      {
        corefilename += scratchDir;
        corefilename += "/";
      }


    corefilename +=  "HYBRIDV2_core_"
        + std::to_string((now->tm_year + 1900)) + "-"
        + std::to_string((now->tm_mon + 1))   + "-"
        + std::to_string(now->tm_mday) + "_" + hostname + "_pid_" + std::to_string(::getpid());

    WriteCoreDump(corefilename.c_str());
    return corefilename;
  }
#else
  std::string dump() {
    return {};
  }
#endif


  void
  warn(
      std::string msg,
      const std::string& file,
      const int line)
  {
#pragma omp critical  (COUT)
    {
      int mpi_is_initialized;
//      MPI_Initialized(&mpi_is_initialized);

      int rank;
//      MPI_Comm_rank(MPI_COMM_WORLD,&rank);

      msg = "\nTEST FAILED ON MPI WORLD RANK " + std::to_string(rank) + " at " + file
          + (
              line >= 0 ?
                  ( std::string(", l.") + std::to_string(line) )
                  : std::string(" ")
          )
          + ": " + msg + "\n";

      fprintf(stderr,msg.c_str());
    }

  }

  void
  terminate(
      std::string msg,
      const char* file,
      const int line)
  {

//    int mpi_is_initialized;
//    MPI_Initialized(&mpi_is_initialized);

//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//    msg = "\nASSERTION FAILED ON MPI WORLD RANK " + std::to_string(rank) + " at " + file
//        + (
//            line >= 0 ?
//                ( std::string(", l.") + std::to_string(line) )
//                : std::string(" ")
//        )
//        + ": " + msg + "\n";

    msg = "\nASSERTION FAILED at " + std::string(file)
        + (
            line >= 0 ?
                ( std::string(", l.") + std::to_string(line) )
                : std::string(" ")
        )
        + ": " + msg + "\n";

    msg += print_stacktrace();

    if (false)//(mpi_is_initialized)
      {
#ifdef DUMP_CORE
        msg += "\n" + std::to_string(rank) + " attempting core dump...\n";
        std::string cf = dump();
        msg += "..." + std::to_string(rank) + " wrote " + cf + ", ";
#endif
        msg += "will attempt mpi shutdown in 30s\n";
        fprintf(stderr,msg.c_str());
        Timeout::sleep_for(std::chrono::seconds(30));
//        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    else
      {
        msg += "\nattempting core dump...\n";
        std::string cf = dump();
        msg += "... wrote " + cf + ", exiting\n";
        fprintf(stderr,msg.c_str());
      }

    exit(1);

  }

} /* namespace MyError */






