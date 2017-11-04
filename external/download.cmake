if(WIN32)
  set(GMP_URL https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/GMP/5.0.1/)
  set(GMP_NAME gmp-all-CGAL-3.9.zip)
  set(GMP_DIR "${CGAL_DIR}/auxiliary/gmp")
  set(MPFR_URL https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/MPFR/3.0.0/)
  set(MPFR_NAME mpfr-all-CGAL-3.9.zip)
  set(MPFR_DIR "${CGAL_DIR}/auxiliary/gmp")

  file(DOWNLOAD "${GMP_URL}/${GMP_NAME}" ${GMP_DIR}/${GMP_NAME})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${GMP_DIR}/${GMP_NAME} WORKING_DIRECTORY ${GMP_DIR})

  file(DOWNLOAD "${MPFR_URL}/${MPFR_NAME}" ${MPFR_DIR}/${MPFR_NAME})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${MPFR_DIR}/${MPFR_NAME} WORKING_DIRECTORY ${MPFR_DIR})
endif()
