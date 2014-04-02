#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GEANT5)


#---Set pathes where to put the libraries, executables and headers------------------------------
set(OUTDIR $ENV{GEANT_OUTPUT_DIR})
if(NOT OUTDIR)
  set(OUTDIR ${PROJECT_SOURCE_DIR})
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTDIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTDIR}/bin)

#----------------------------------------------------------------------------
# Find ROOT and setup the libraries
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIR})
include_directories(/home/ramanNew/BOOST/nt2/nt2-install/include)
#----------------------------------------------------------------------------
if(CUSTOM_ALLOCATOR)
  find_package(CAllocator)
endif()

#----------------------------------------------------------------------------
# Add subdirectories
#
add_custom_target(_headers ALL)
add_subdirectory (geom_vec_tests)

