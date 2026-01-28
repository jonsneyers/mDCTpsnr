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

/// Includes
#include "thread.hpp"
#include "std/stdarg.hpp"
#include "std/assert.hpp"
#include "std/string.hpp"
///

/// Thread::Thread
Thread::Thread(void)
  : m_pThreads(NULL)
{ }
///

/// Thread::~Thread
Thread::~Thread(void)
{
  int i;
  
  for(i = 0;i < m_iTotal;i++) {
    if (m_pThreads[i].m_bRunning) {
      m_pThreads[i].m_bDie = true;
      sem_post(&m_pThreads[i].m_MayRun);
      pthread_join(m_pThreads[i].m_Thread,NULL);
    }
  }
  delete[] m_pThreads;
}
///

/// Thread::EntryPoint
void *Thread::EntryPoint(void *arg)
{
  struct ThreadStarter *ts = (struct ThreadStarter *)arg;

  do {
    sem_wait(&ts->m_MayRun);
    if (ts->m_bDie)
      break;
    try {
      ts->m_That->Run(ts->m_iCPU,ts->m_iCores);
    } catch(...) {
      // Ignore. Yuck!
    }
    sem_post(&ts->m_Done);
  } while(true);

  return ts;
}
///

/// Thread::SplitWork
// Launch the side-threads for the given number of cores, run the own computation
// on this core, return when all are done.
void Thread::SplitWork(int cores)
{
  int i;

  if (m_pThreads == NULL) {
    m_pThreads = new struct ThreadStarter[cores];
    m_iTotal   = cores;
    for(i = 0;i < cores;i++) {
      m_pThreads[i].m_That     = this;
      m_pThreads[i].m_iCPU     = i;
      m_pThreads[i].m_iCores   = cores;
      m_pThreads[i].m_bRunning = false;
      m_pThreads[i].m_bDie     = false;
      sem_init(&m_pThreads[i].m_MayRun,0,0);
      sem_init(&m_pThreads[i].m_Done  ,0,0);
      if (i != 0) {
	// The first thread is run on this core.
	if (pthread_create(&m_pThreads[i].m_Thread,NULL,&Thread::EntryPoint,&m_pThreads[i]) == 0) {
	  m_pThreads[i].m_bRunning = true;
	}
      }
    }
  } else {
    assert(cores == m_iTotal);
  }

  // Launch the threads.
  for(i = 0;i < cores;i++) {
    if (m_pThreads[i].m_bRunning)
      sem_post(&m_pThreads[i].m_MayRun);
  }

  // Do all non-running threads myself.
  for(i = 0;i < cores;i++) {
    if (!m_pThreads[i].m_bRunning) {
      Run(i,m_iTotal);
    }
  }

  // Now collect the threads again.
  for(i = 0;i < cores;i++) {
    if (m_pThreads[i].m_bRunning) {
      sem_wait(&m_pThreads[i].m_Done);
    }
  }
}
///
