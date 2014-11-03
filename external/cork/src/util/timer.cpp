// +-------------------------------------------------------------------------
// | timer.cpp
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#include "prelude.h"

#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>
void gettimeofday(struct timeval* t,void* timezone)
{       struct _timeb timebuffer;
        _ftime( &timebuffer );
        t->tv_sec=timebuffer.time;
        t->tv_usec=1000*timebuffer.millitm;
}
#else
#include <sys/time.h>
#endif





Timer::Timer()
{
    start();
}
Timer::~Timer()
{
}

inline
double timeDiff(const timeval &t0, const timeval &t1)
{
    double ms = (t1.tv_sec - t0.tv_sec) * 1000.0
              + (t1.tv_usec - t0.tv_usec) / 1000.0;
    return ms;
}

void Timer::start()
{
    gettimeofday(&init, NULL);
    prev = init;
    last_lap_duration = ellapsed_time_duration = 0.0;
    isRunning = true;
}
// returns lap time in milliseconds
double Timer::lap()
{
    if(!isRunning)
        return 0.0;
    
    timeval now;
    gettimeofday(&now, NULL);
    
    last_lap_duration       = timeDiff(prev, now);
    return last_lap_duration;
}
// returns total ellapsed time in milliseconds
double Timer::stop()
{
    if(!isRunning)
        return 0.0;
    
    timeval now;
    gettimeofday(&now, NULL);
    
    last_lap_duration       = timeDiff(prev, now);
    ellapsed_time_duration  = timeDiff(init, now);
    
    isRunning = false;
    return ellapsed_time_duration;
}
// re-retrieve values
double Timer::lastLap() const
{
    return last_lap_duration;
}
double Timer::ellapsed() const
{
    return ellapsed_time_duration;
}
