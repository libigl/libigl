// +-------------------------------------------------------------------------
// | log.cpp
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

#include <cstdlib>
using std::atexit;

#include <ctime>
#include <fstream>
using std::ofstream;
using std::endl;

namespace {

ofstream error_log_stream;

void on_exit()
{
    error_log_stream << "Ending error logging at " << endl;
    std::time_t time_var = std::time(NULL);
    error_log_stream << std::ctime(&time_var) << endl;
    error_log_stream.close();
}

}

void logInit()
{
    error_log_stream.open("error_log.txt", std::ios_base::app | std::ios_base::out );
    error_log_stream << "Begining error logging at " << endl;
    std::time_t time_var = std::time(NULL);
    error_log_stream << std::ctime(&time_var);// << endl;
    
    atexit(on_exit);
}

// exposes the error log for writing
std::ostream& err()
{
    static bool first = true;
    if(first) {
        logInit();
        first = false;
    }
    return error_log_stream;
}


