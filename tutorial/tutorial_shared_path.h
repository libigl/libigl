#ifndef tutorial_shared_path_h_included
#define tutorial_shared_path_h_included

// convert the argument to a string constant
#define STRINGIZE2(s) #s
// need an extra level of macro indirection for preprocessor token stringification to work
#define STRINGIZE(s) STRINGIZE2(s)

#include <string>

#ifndef TUTORIAL_SHARED_PATH
const std::string tutorial_shared_path("../shared");
#else
const std::string tutorial_shared_path( STRINGIZE(TUTORIAL_SHARED_PATH) );
#endif

#endif
