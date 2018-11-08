#include <gtest/gtest.h>

#ifndef NDEBUG
#ifdef __linux__
#include <fenv.h>
#endif
#endif


int main(int argc, char **argv) {
    
#ifndef NDEBUG
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
#endif     
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
