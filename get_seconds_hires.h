#ifndef IGL_GET_SECONDS_HIRES_H
#define IGL_GET_SECONDS_HIRES_H

namespace igl
{
	// Return the current time in seconds using performance counters
	inline double get_seconds_hires();
}

//Implementation
#if _WIN32
#  include <windows.h>
inline double igl::get_seconds_hires()
{
	LARGE_INTEGER li_freq, li_current;
	const bool ret = QueryPerformanceFrequency(&li_freq);
	const bool ret2 = QueryPerformanceCounter(&li_current);
	assert(ret && ret2);
	assert(li_freq.QuadPart > 0);
	return double(li_current.QuadPart) / double(li_freq.QuadPart);
}
#else
#  include "get_seconds.h"
inline double igl::get_seconds_hires()
{
	// Sorry I've no idea how performance counters work on Mac...
	return igl::get_seconds();
}
#endif
#endif
