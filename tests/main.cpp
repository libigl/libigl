////////////////////////////////////////////////////////////////////////////////
// Keep this file empty, and implement unit tests in separate compilation units!
////////////////////////////////////////////////////////////////////////////////

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>


// TODO: Fix floating point exceptions raised in debug mode before re-enabling this.
// #ifndef NDEBUG
// #ifdef __linux__
// #include <fenv.h>
// #endif
// #endif

// #ifndef NDEBUG
// #ifdef __linux__
// void beforeMain (void) __attribute__((constructor));
// void beforeMain (void)
// {
//     feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
// }
// #endif
// #endif
