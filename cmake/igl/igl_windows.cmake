if(MSVC)
	# https://github.com/mozilla/sccache/issues/242
	if(CMAKE_CXX_COMPILER_LAUNCHER STREQUAL "sccache")
		string(REGEX REPLACE "/Z[iI7]" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Z7")
	endif()
endif()
