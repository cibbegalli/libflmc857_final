#libFL CMAKE

cmake_minimum_required(VERSION 3.1)

project(libFL)

set(CMAKE_VERBOSE_MAKEFILE "ON")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Threads REQUIRED)

message(STATUS ${CMAKE_MODULE_PATH})

find_package( CBLAS REQUIRED )
find_package( LAPACK REQUIRED )
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

option (USE_LAPACKE "Use LAPACKE C interface to LAPACK" ON)
option (USE_LAPACKE "Use BLAS C interface to BLAS" ON)
option (USE_OPENMP "Use OpenMP for multithreading" ON)
 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
 
include_directories(include src .)
message(Compiler: "${CMAKE_CXX_COMPILER_ID}")

message(Compiler: "${CMAKE_C_COMPILER_ID}")
 
aux_source_directory(${libFL_SOURCE_DIR}/src LIBFL_SRC)
 
file(GLOB LIBFL_INCLUDE
    ${libfl_SOURCE_DIR}/include/*.h)
 

message(STATUS "LIBFL library")
foreach( libflsourcefile ${LIBFL_SRC} )
    message(STATUS "\t${libflsourcefile}")
endforeach(libflsourcefile ${LIBFL_SRC})
 
add_library(libfl STATIC ${LIBFL_SRC} ${LIBFL_INCLUDE})

file(GLOB LIBFL_DEMO

    demo/*.cpp
    )

install(FILES ${LIBFL_INCLUDE}
        DESTINATION include
        )
 
message(STATUS "Demo files")
foreach( demosourcefile ${LIBFL_DEMO} )
    get_filename_component(demofilename ${demosourcefile} NAME)
    string( REPLACE ".cpp" "" demoname ${demofilename} )

    message(STATUS "\t${demosourcefile}")

    add_executable( ${demoname} ${demosourcefile})
    target_link_libraries( ${demoname} libfl ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES} stdc++)

    #install(TARGETS ${demoname} RUNTIME DESTINATION bin)

    #remove demo from make all
    #set_target_properties(${demoname} PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
endforeach( demosourcefile ${LIBFL_DEMO} )

install (TARGETS libfl
         ARCHIVE DESTINATION libfl
         LIBRARY DESTINATION libfl
        )#RUNTIME DESTINATION bin)
