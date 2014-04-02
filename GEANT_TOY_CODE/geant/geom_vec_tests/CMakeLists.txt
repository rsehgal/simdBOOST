set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})

include_directories(inc)
include_directories(/home/ramanNew/BOOST/nt2/nt2-install/include)

find_package(TBB)
include_directories(${TBB_INCLUDE_DIR})

file(GLOB sources src/TGeoBBox_v.cxx src/TGeoTube_v.cxx )
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)


if( USE_BOOST_SIMD )
  add_definitions( "-DBOOST_SIMD" )
  #set( BASEVCFLAGS "-std=c++11 -fabi-version=6" )
endif()

if( USEAVX )
  	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} -mavx")
else()
   	set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS} -msse4.1")
endif()


add_library(VecGeom SHARED ${sources} ${headers})
add_library(VecGeom_static STATIC ${sources} ${headers})
target_link_libraries(VecGeom ${ROOT_LIBRARIES} -lEG -lGeom)

add_custom_target(VectorizedGeom DEPENDS VecGeom VecGeom_static)


#----------------------------------------------------------------------------
# Add some benchmarks and standard tests if wanted

#file(GLOB benchmarksources src/tests/TestSIMD.cpp)
file(GLOB benchmarksources src/BenchmarkFramework/BenchTGeoBBox_v.cxx
			   src/BenchmarkFramework/BenchTGeoTube_v.cxx)
set_source_files_properties(${benchmarksources} PROPERTIES COMPILE_FLAGS "${BASEFLAGS}")

    foreach( benchmarksrc  ${benchmarksources} )
        get_filename_component( BENCHEXEC  ${benchmarksrc} NAME_WE)
	message(STATUS "building benchmark ${BENCHEXEC}")
	add_executable( ${BENCHEXEC} ${benchmarksrc} ${headers} )
	target_link_libraries( ${BENCHEXEC} ${ROOT_LIBRARIES} -lGeom -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lVecGeom_static ${TBB_LIBRARIES})
	add_dependencies( ${BENCHEXEC} VectorizedGeom )
    endforeach( benchmarksrc )



#----------------------------------------------------------------------------
 
