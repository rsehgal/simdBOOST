# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ramanNew/GEANT_TOY_CODE/geant

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ramanNew/GEANT_TOY_CODE/geant

# Include any dependencies generated for this target.
include geom_vec_tests/CMakeFiles/VecGeom.dir/depend.make

# Include the progress variables for this target.
include geom_vec_tests/CMakeFiles/VecGeom.dir/progress.make

# Include the compile flags for this target's objects.
include geom_vec_tests/CMakeFiles/VecGeom.dir/flags.make

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o: geom_vec_tests/CMakeFiles/VecGeom.dir/flags.make
geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o: geom_vec_tests/src/TGeoBBox_v.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ramanNew/GEANT_TOY_CODE/geant/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/gcc   $(CXX_DEFINES) $(CXX_FLAGS)  -mavx -o CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o -c /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/TGeoBBox_v.cxx

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.i"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/gcc  $(CXX_DEFINES) $(CXX_FLAGS)  -mavx -E /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/TGeoBBox_v.cxx > CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.i

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.s"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/gcc  $(CXX_DEFINES) $(CXX_FLAGS)  -mavx -S /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/TGeoBBox_v.cxx -o CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.s

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.requires:
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.requires

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.provides: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.requires
	$(MAKE) -f geom_vec_tests/CMakeFiles/VecGeom.dir/build.make geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.provides.build
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.provides

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.provides.build: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o: geom_vec_tests/CMakeFiles/VecGeom.dir/flags.make
geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o: geom_vec_tests/src/TGeoTube_v.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ramanNew/GEANT_TOY_CODE/geant/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/gcc   $(CXX_DEFINES) $(CXX_FLAGS)  -mavx -o CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o -c /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/TGeoTube_v.cxx

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.i"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/gcc  $(CXX_DEFINES) $(CXX_FLAGS)  -mavx -E /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/TGeoTube_v.cxx > CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.i

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.s"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/gcc  $(CXX_DEFINES) $(CXX_FLAGS)  -mavx -S /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/TGeoTube_v.cxx -o CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.s

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.requires:
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.requires

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.provides: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.requires
	$(MAKE) -f geom_vec_tests/CMakeFiles/VecGeom.dir/build.make geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.provides.build
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.provides

geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.provides.build: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o

# Object files for target VecGeom
VecGeom_OBJECTS = \
"CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o" \
"CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o"

# External object files for target VecGeom
VecGeom_EXTERNAL_OBJECTS =

lib/libVecGeom.so: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o
lib/libVecGeom.so: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o
lib/libVecGeom.so: geom_vec_tests/CMakeFiles/VecGeom.dir/build.make
lib/libVecGeom.so: geom_vec_tests/CMakeFiles/VecGeom.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../lib/libVecGeom.so"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/VecGeom.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
geom_vec_tests/CMakeFiles/VecGeom.dir/build: lib/libVecGeom.so
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/build

geom_vec_tests/CMakeFiles/VecGeom.dir/requires: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoBBox_v.cxx.o.requires
geom_vec_tests/CMakeFiles/VecGeom.dir/requires: geom_vec_tests/CMakeFiles/VecGeom.dir/src/TGeoTube_v.cxx.o.requires
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/requires

geom_vec_tests/CMakeFiles/VecGeom.dir/clean:
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && $(CMAKE_COMMAND) -P CMakeFiles/VecGeom.dir/cmake_clean.cmake
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/clean

geom_vec_tests/CMakeFiles/VecGeom.dir/depend:
	cd /home/ramanNew/GEANT_TOY_CODE/geant && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ramanNew/GEANT_TOY_CODE/geant /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests /home/ramanNew/GEANT_TOY_CODE/geant /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/CMakeFiles/VecGeom.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : geom_vec_tests/CMakeFiles/VecGeom.dir/depend

