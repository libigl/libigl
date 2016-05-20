// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PROGRESS_H
#define IGL_PROGRESS_H
#include "igl_inline.h"

namespace igl
{
  // Class defining interface for support interruptable task.
  //
  // The default class is a no-op, relying on the class user to give it
  // application specific implementations.  Typically, an application has a
  // singleton progress class that is then passed onto various functions which
  // use it to check if they should be interrupted.
  class Progress
  {
    // Default constructor
    Progress() {}

    // Signal the start of an interruptible task with an optional name.
    // In general, this interface supports nested tasks. A matching call to
    // endTask() is required, unless the progress was interrupted.
    // If startTask() returns false, then it means that the user already
    // interrupted and we should not continue, but endTask() still needs to be
    // called.
    bool startTask(const char* = nullptr) { return true; }

    // Signal the end of an interruptible task. Alwatys call this if you called
    // startTask() (regardless of its return value).
    void endTask() {}

    // Update the progress and check if the task should be aborted.
    //
    // Inputs:
    //	 percent  an optional (when >= 0) percentage indicating
    //		  the fraction of the task that has been completed
    //
    // Returns true if task should be interrupted.
    IGL_INLINE bool wasInterrupted(int = -1)
      { return false; }

  };

  typedef Progress NullProgress;

  // This method allows Progress::wasInterrupted to be compiled out when
  // client code only has a pointer (vs reference) to the progress class.
  //
  // NOTE: This is a free-standing function since C++ doesn't allow for
  // partial template specialization (in client code of the progress class).
  template <typename T>
  IGL_INLINE bool was_interrupted(T* i, int percent = -1)
    { return i && i->wasInterrupted(percent); }

  // Specialization for NullProgress to ensure it gets compiled out
  template<>
  IGL_INLINE bool was_interrupted<NullProgress>(NullProgress*, int)
    { return false; }

  // Class defining the scope of an interruptable task
  template <typename ProgressT>
  class ProgressTask
  {
  public:
    // Start an interruptable task with optional name
    ProgressTask(ProgressT* progress, const char* name = nullptr)
      : myProgress(progress)
      , myIsRunning(progress ? progress->startTask(name) : false)
    {
    }

    ~ProgressTask()
    {
      if (myProgress)
	myProgress->endTask();
    }

    // Update the progress and check if the task should be aborted.
    //
    // Inputs:
    //	 percent  an optional (when >= 0) percentage indicating
    //		  the fraction of the task that has been completed
    //
    // Returns true if task should be interrupted.
    IGL_INLINE bool wasInterrupted(int percent = -1)
    {
	if (myIsRunning && igl::was_interrupted(myProgress, percent))
	    myIsRunning = false;
	return !myIsRunning;
    }

  private:
    ProgressT* myProgress;
    bool myIsRunning;
  };

  // Specialization for NullProgress to ensure it gets compiled out
  template <>
  class ProgressTask<NullProgress>
  {
  public:
    // Start an interruptable task with optional name
    ProgressTask(NullProgress*, const char* name) {}

    ~ProgressTask() {}

    IGL_INLINE bool wasInterrupted(int = -1) { return false; }
  };

}
#endif // IGL_PROGRESS_H
