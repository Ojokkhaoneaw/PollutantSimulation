# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.3.4\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.3.4\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\JirayuH\PollutantSimulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\JirayuH\PollutantSimulation\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/compute_poisson.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/compute_poisson.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/compute_poisson.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/compute_poisson.dir/flags.make

CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj: CMakeFiles/compute_poisson.dir/flags.make
CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj: ../compute_poisson.cpp
CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj: CMakeFiles/compute_poisson.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\JirayuH\PollutantSimulation\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.4\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj -MF CMakeFiles\compute_poisson.dir\compute_poisson.cpp.obj.d -o CMakeFiles\compute_poisson.dir\compute_poisson.cpp.obj -c C:\Users\JirayuH\PollutantSimulation\compute_poisson.cpp

CMakeFiles/compute_poisson.dir/compute_poisson.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compute_poisson.dir/compute_poisson.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.4\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\JirayuH\PollutantSimulation\compute_poisson.cpp > CMakeFiles\compute_poisson.dir\compute_poisson.cpp.i

CMakeFiles/compute_poisson.dir/compute_poisson.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compute_poisson.dir/compute_poisson.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.4\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\JirayuH\PollutantSimulation\compute_poisson.cpp -o CMakeFiles\compute_poisson.dir\compute_poisson.cpp.s

# Object files for target compute_poisson
compute_poisson_OBJECTS = \
"CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj"

# External object files for target compute_poisson
compute_poisson_EXTERNAL_OBJECTS =

libcompute_poisson.a: CMakeFiles/compute_poisson.dir/compute_poisson.cpp.obj
libcompute_poisson.a: CMakeFiles/compute_poisson.dir/build.make
libcompute_poisson.a: CMakeFiles/compute_poisson.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\JirayuH\PollutantSimulation\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcompute_poisson.a"
	$(CMAKE_COMMAND) -P CMakeFiles\compute_poisson.dir\cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\compute_poisson.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/compute_poisson.dir/build: libcompute_poisson.a
.PHONY : CMakeFiles/compute_poisson.dir/build

CMakeFiles/compute_poisson.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\compute_poisson.dir\cmake_clean.cmake
.PHONY : CMakeFiles/compute_poisson.dir/clean

CMakeFiles/compute_poisson.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\JirayuH\PollutantSimulation C:\Users\JirayuH\PollutantSimulation C:\Users\JirayuH\PollutantSimulation\cmake-build-debug C:\Users\JirayuH\PollutantSimulation\cmake-build-debug C:\Users\JirayuH\PollutantSimulation\cmake-build-debug\CMakeFiles\compute_poisson.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/compute_poisson.dir/depend

