# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build"

# Include any dependencies generated for this target.
include CMakeFiles/confMap.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/confMap.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/confMap.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/confMap.dir/flags.make

CMakeFiles/confMap.dir/confMap.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/confMap.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/confMap.cpp
CMakeFiles/confMap.dir/confMap.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/confMap.dir/confMap.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/confMap.cpp.o -MF CMakeFiles/confMap.dir/confMap.cpp.o.d -o CMakeFiles/confMap.dir/confMap.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/confMap.cpp"

CMakeFiles/confMap.dir/confMap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/confMap.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/confMap.cpp" > CMakeFiles/confMap.dir/confMap.cpp.i

CMakeFiles/confMap.dir/confMap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/confMap.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/confMap.cpp" -o CMakeFiles/confMap.dir/confMap.cpp.s

CMakeFiles/confMap.dir/Core/Edge.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Edge.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Edge.cpp
CMakeFiles/confMap.dir/Core/Edge.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/confMap.dir/Core/Edge.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Edge.cpp.o -MF CMakeFiles/confMap.dir/Core/Edge.cpp.o.d -o CMakeFiles/confMap.dir/Core/Edge.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Edge.cpp"

CMakeFiles/confMap.dir/Core/Edge.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Edge.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Edge.cpp" > CMakeFiles/confMap.dir/Core/Edge.cpp.i

CMakeFiles/confMap.dir/Core/Edge.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Edge.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Edge.cpp" -o CMakeFiles/confMap.dir/Core/Edge.cpp.s

CMakeFiles/confMap.dir/Core/Face.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Face.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Face.cpp
CMakeFiles/confMap.dir/Core/Face.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/confMap.dir/Core/Face.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Face.cpp.o -MF CMakeFiles/confMap.dir/Core/Face.cpp.o.d -o CMakeFiles/confMap.dir/Core/Face.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Face.cpp"

CMakeFiles/confMap.dir/Core/Face.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Face.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Face.cpp" > CMakeFiles/confMap.dir/Core/Face.cpp.i

CMakeFiles/confMap.dir/Core/Face.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Face.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Face.cpp" -o CMakeFiles/confMap.dir/Core/Face.cpp.s

CMakeFiles/confMap.dir/Core/FException.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/FException.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/FException.cpp
CMakeFiles/confMap.dir/Core/FException.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/confMap.dir/Core/FException.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/FException.cpp.o -MF CMakeFiles/confMap.dir/Core/FException.cpp.o.d -o CMakeFiles/confMap.dir/Core/FException.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/FException.cpp"

CMakeFiles/confMap.dir/Core/FException.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/FException.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/FException.cpp" > CMakeFiles/confMap.dir/Core/FException.cpp.i

CMakeFiles/confMap.dir/Core/FException.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/FException.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/FException.cpp" -o CMakeFiles/confMap.dir/Core/FException.cpp.s

CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/HalfEdge.cpp
CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o -MF CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o.d -o CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/HalfEdge.cpp"

CMakeFiles/confMap.dir/Core/HalfEdge.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/HalfEdge.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/HalfEdge.cpp" > CMakeFiles/confMap.dir/Core/HalfEdge.cpp.i

CMakeFiles/confMap.dir/Core/HalfEdge.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/HalfEdge.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/HalfEdge.cpp" -o CMakeFiles/confMap.dir/Core/HalfEdge.cpp.s

CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/OBJFileReader.cpp
CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o -MF CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o.d -o CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/OBJFileReader.cpp"

CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/OBJFileReader.cpp" > CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.i

CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/OBJFileReader.cpp" -o CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.s

CMakeFiles/confMap.dir/Core/Parser.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Parser.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Parser.cpp
CMakeFiles/confMap.dir/Core/Parser.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/confMap.dir/Core/Parser.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Parser.cpp.o -MF CMakeFiles/confMap.dir/Core/Parser.cpp.o.d -o CMakeFiles/confMap.dir/Core/Parser.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Parser.cpp"

CMakeFiles/confMap.dir/Core/Parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Parser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Parser.cpp" > CMakeFiles/confMap.dir/Core/Parser.cpp.i

CMakeFiles/confMap.dir/Core/Parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Parser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Parser.cpp" -o CMakeFiles/confMap.dir/Core/Parser.cpp.s

CMakeFiles/confMap.dir/Core/Point.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Point.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Point.cpp
CMakeFiles/confMap.dir/Core/Point.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/confMap.dir/Core/Point.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Point.cpp.o -MF CMakeFiles/confMap.dir/Core/Point.cpp.o.d -o CMakeFiles/confMap.dir/Core/Point.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Point.cpp"

CMakeFiles/confMap.dir/Core/Point.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Point.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Point.cpp" > CMakeFiles/confMap.dir/Core/Point.cpp.i

CMakeFiles/confMap.dir/Core/Point.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Point.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Point.cpp" -o CMakeFiles/confMap.dir/Core/Point.cpp.s

CMakeFiles/confMap.dir/Core/Solid.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Solid.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Solid.cpp
CMakeFiles/confMap.dir/Core/Solid.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/confMap.dir/Core/Solid.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Solid.cpp.o -MF CMakeFiles/confMap.dir/Core/Solid.cpp.o.d -o CMakeFiles/confMap.dir/Core/Solid.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Solid.cpp"

CMakeFiles/confMap.dir/Core/Solid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Solid.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Solid.cpp" > CMakeFiles/confMap.dir/Core/Solid.cpp.i

CMakeFiles/confMap.dir/Core/Solid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Solid.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Solid.cpp" -o CMakeFiles/confMap.dir/Core/Solid.cpp.s

CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/SolidDelegate.cpp
CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o -MF CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o.d -o CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/SolidDelegate.cpp"

CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/SolidDelegate.cpp" > CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.i

CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/SolidDelegate.cpp" -o CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.s

CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/StringTokenizer.cpp
CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o -MF CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o.d -o CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/StringTokenizer.cpp"

CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/StringTokenizer.cpp" > CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.i

CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/StringTokenizer.cpp" -o CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.s

CMakeFiles/confMap.dir/Core/TopologyException.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/TopologyException.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/TopologyException.cpp
CMakeFiles/confMap.dir/Core/TopologyException.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/confMap.dir/Core/TopologyException.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/TopologyException.cpp.o -MF CMakeFiles/confMap.dir/Core/TopologyException.cpp.o.d -o CMakeFiles/confMap.dir/Core/TopologyException.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/TopologyException.cpp"

CMakeFiles/confMap.dir/Core/TopologyException.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/TopologyException.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/TopologyException.cpp" > CMakeFiles/confMap.dir/Core/TopologyException.cpp.i

CMakeFiles/confMap.dir/Core/TopologyException.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/TopologyException.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/TopologyException.cpp" -o CMakeFiles/confMap.dir/Core/TopologyException.cpp.s

CMakeFiles/confMap.dir/Core/Trait.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Trait.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Trait.cpp
CMakeFiles/confMap.dir/Core/Trait.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/confMap.dir/Core/Trait.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Trait.cpp.o -MF CMakeFiles/confMap.dir/Core/Trait.cpp.o.d -o CMakeFiles/confMap.dir/Core/Trait.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Trait.cpp"

CMakeFiles/confMap.dir/Core/Trait.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Trait.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Trait.cpp" > CMakeFiles/confMap.dir/Core/Trait.cpp.i

CMakeFiles/confMap.dir/Core/Trait.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Trait.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Trait.cpp" -o CMakeFiles/confMap.dir/Core/Trait.cpp.s

CMakeFiles/confMap.dir/Core/Vertex.cpp.o: CMakeFiles/confMap.dir/flags.make
CMakeFiles/confMap.dir/Core/Vertex.cpp.o: /home/manjaro/Documents/Advanced\ Computer\ Graphics\ 1/Project\ 3/Spherical-Conformal-Mapping/Test/Core/Vertex.cpp
CMakeFiles/confMap.dir/Core/Vertex.cpp.o: CMakeFiles/confMap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/confMap.dir/Core/Vertex.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/confMap.dir/Core/Vertex.cpp.o -MF CMakeFiles/confMap.dir/Core/Vertex.cpp.o.d -o CMakeFiles/confMap.dir/Core/Vertex.cpp.o -c "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Vertex.cpp"

CMakeFiles/confMap.dir/Core/Vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/confMap.dir/Core/Vertex.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Vertex.cpp" > CMakeFiles/confMap.dir/Core/Vertex.cpp.i

CMakeFiles/confMap.dir/Core/Vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/confMap.dir/Core/Vertex.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/Core/Vertex.cpp" -o CMakeFiles/confMap.dir/Core/Vertex.cpp.s

# Object files for target confMap
confMap_OBJECTS = \
"CMakeFiles/confMap.dir/confMap.cpp.o" \
"CMakeFiles/confMap.dir/Core/Edge.cpp.o" \
"CMakeFiles/confMap.dir/Core/Face.cpp.o" \
"CMakeFiles/confMap.dir/Core/FException.cpp.o" \
"CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o" \
"CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o" \
"CMakeFiles/confMap.dir/Core/Parser.cpp.o" \
"CMakeFiles/confMap.dir/Core/Point.cpp.o" \
"CMakeFiles/confMap.dir/Core/Solid.cpp.o" \
"CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o" \
"CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o" \
"CMakeFiles/confMap.dir/Core/TopologyException.cpp.o" \
"CMakeFiles/confMap.dir/Core/Trait.cpp.o" \
"CMakeFiles/confMap.dir/Core/Vertex.cpp.o"

# External object files for target confMap
confMap_EXTERNAL_OBJECTS =

confMap: CMakeFiles/confMap.dir/confMap.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Edge.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Face.cpp.o
confMap: CMakeFiles/confMap.dir/Core/FException.cpp.o
confMap: CMakeFiles/confMap.dir/Core/HalfEdge.cpp.o
confMap: CMakeFiles/confMap.dir/Core/OBJFileReader.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Parser.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Point.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Solid.cpp.o
confMap: CMakeFiles/confMap.dir/Core/SolidDelegate.cpp.o
confMap: CMakeFiles/confMap.dir/Core/StringTokenizer.cpp.o
confMap: CMakeFiles/confMap.dir/Core/TopologyException.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Trait.cpp.o
confMap: CMakeFiles/confMap.dir/Core/Vertex.cpp.o
confMap: CMakeFiles/confMap.dir/build.make
confMap: CMakeFiles/confMap.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_15) "Linking CXX executable confMap"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/confMap.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/confMap.dir/build: confMap
.PHONY : CMakeFiles/confMap.dir/build

CMakeFiles/confMap.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/confMap.dir/cmake_clean.cmake
.PHONY : CMakeFiles/confMap.dir/clean

CMakeFiles/confMap.dir/depend:
	cd "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test" "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test" "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build" "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build" "/home/manjaro/Documents/Advanced Computer Graphics 1/Project 3/Spherical-Conformal-Mapping/Test/build/CMakeFiles/confMap.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : CMakeFiles/confMap.dir/depend

