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
include geom_vec_tests/CMakeFiles/TestSIMD.dir/depend.make

# Include the progress variables for this target.
include geom_vec_tests/CMakeFiles/TestSIMD.dir/progress.make

# Include the compile flags for this target's objects.
include geom_vec_tests/CMakeFiles/TestSIMD.dir/flags.make

geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o: geom_vec_tests/CMakeFiles/TestSIMD.dir/flags.make
geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o: geom_vec_tests/src/tests/TestSIMD.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ramanNew/GEANT_TOY_CODE/geant/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o -c /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/tests/TestSIMD.cpp

geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.i"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/tests/TestSIMD.cpp > CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.i

geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.s"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && /home/ramanNew/GCC-4.7.2/gcc-4.8.1/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/src/tests/TestSIMD.cpp -o CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.s

geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.requires:
.PHONY : geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.requires

geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.provides: geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.requires
	$(MAKE) -f geom_vec_tests/CMakeFiles/TestSIMD.dir/build.make geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.provides.build
.PHONY : geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.provides

geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.provides.build: geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o

# Object files for target TestSIMD
TestSIMD_OBJECTS = \
"CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o"

# External object files for target TestSIMD
TestSIMD_EXTERNAL_OBJECTS =

bin/TestSIMD: geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o
bin/TestSIMD: geom_vec_tests/CMakeFiles/TestSIMD.dir/build.make
bin/TestSIMD: /home/ramanNew/GEANT_PREREQ/tbb/lib/libtbb.so
bin/TestSIMD: geom_vec_tests/CMakeFiles/TestSIMD.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/TestSIMD"
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestSIMD.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
geom_vec_tests/CMakeFiles/TestSIMD.dir/build: bin/TestSIMD
.PHONY : geom_vec_tests/CMakeFiles/TestSIMD.dir/build

geom_vec_tests/CMakeFiles/TestSIMD.dir/requires: geom_vec_tests/CMakeFiles/TestSIMD.dir/src/tests/TestSIMD.cpp.o.requires
.PHONY : geom_vec_tests/CMakeFiles/TestSIMD.dir/requires

geom_vec_tests/CMakeFiles/TestSIMD.dir/clean:
	cd /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests && $(CMAKE_COMMAND) -P CMakeFiles/TestSIMD.dir/cmake_clean.cmake
.PHONY : geom_vec_tests/CMakeFiles/TestSIMD.dir/clean

geom_vec_tests/CMakeFiles/TestSIMD.dir/depend:
	cd /home/ramanNew/GEANT_TOY_CODE/geant && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ramanNew/GEANT_TOY_CODE/geant /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests /home/ramanNew/GEANT_TOY_CODE/geant /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests /home/ramanNew/GEANT_TOY_CODE/geant/geom_vec_tests/CMakeFiles/TestSIMD.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : geom_vec_tests/CMakeFiles/TestSIMD.dir/depend

