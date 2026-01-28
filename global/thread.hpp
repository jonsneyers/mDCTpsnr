/************************************************************************************
 **  Copyright (C) 2009 University of Stuttgart, Thomas Richter                    **
 **                                                                                **
 **  This software is provided 'as-is', without any express or implied             **
 **  warranty.  In no event will the authors be held liable for any damages        **
 **  arising from the use of this software.                                        **
 **                                                                                **
 **  Permission is granted to anyone to use this software for any purpose,         **
 **  including commercial applications, and to alter it and redistribute it        **
 **  freely, subject to the following restrictions:                                **
 **                                                                                **
 **  1. The origin of this software must not be misrepresented; you must not       **
 **     claim that you wrote the original software. If you use this software       **
 **     in a product, an acknowledgment in the product documentation would be      **
 **     appreciated but is not required.                                           **
 **  2. Altered source versions must be plainly marked as such, and must not be    **
 **     misrepresented as being the original software.                             **
 **  3. This notice may not be removed or altered from any source distribution.    **
 **                                                                                **
 **     Thomas Richter                                                             **
 **     thomas.richter@iis.fraunhofer.de                                           **
 **                                                                                **
 ************************************************************************************/

#ifndef GLOBAL_THREAD_HPP
#define GLOBAL_THREAD_HPP

/// Includes
#include "global/types.hpp"
#include "std/stdarg.hpp"
#include <pthread.h>
#include <semaphore.h>
///

/// class Thread
// This class represents the base class for an activity that can be threaded.
class Thread {
  //
  // The entry point to the thread used by pthread_create.
  static void *EntryPoint(void *arg);
  //
  // This helper function is presented as argument to the thread.
  struct ThreadStarter {
    // 
    // Pointer to the thread class to which this thread belongs.
    class Thread *m_That;
    //
    // The argument to be passed to the "real" entry point, namely the
    // CPU # that is started here.
    int           m_iCPU;
    //
    // The number of cores to run this on.
    int           m_iCores;
    //
    // True in case this thread is running.
    bool          m_bRunning;
    //
    // True in case the thread must die.
    bool          m_bDie;
    //
    // The thread ID
    pthread_t     m_Thread;
    //
    // This semaphore indicates whether the thread may run.
    sem_t         m_MayRun;
    //
    // This semaphore indicates whether the job completed.
    sem_t         m_Done;
    //
  }              *m_pThreads;
  //
  // Number of cores.
  int             m_iTotal;
  //
protected:
  //
  // This is the user function that is called to run the side-thread.
  // Arguments are the CPU# that is run, the total number of CPUs, and
  // additional args.
  virtual void Run(int cpu,int cores) = 0;
  //
public:
  Thread(void);
  virtual ~Thread(void);
  //
  // Launch the side-threads for the given number of cores, run the own computation
  // on this core, return when all are done.
  void SplitWork(int cores);
  //
};
///

///
#endif
