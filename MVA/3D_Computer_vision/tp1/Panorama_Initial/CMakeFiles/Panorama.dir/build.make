# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_SOURCE_DIR = /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial

# Include any dependencies generated for this target.
include CMakeFiles/Panorama.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Panorama.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Panorama.dir/flags.make

CMakeFiles/Panorama.dir/Panorama.cpp.o: CMakeFiles/Panorama.dir/flags.make
CMakeFiles/Panorama.dir/Panorama.cpp.o: Panorama.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Panorama.dir/Panorama.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Panorama.dir/Panorama.cpp.o -c /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial/Panorama.cpp

CMakeFiles/Panorama.dir/Panorama.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Panorama.dir/Panorama.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial/Panorama.cpp > CMakeFiles/Panorama.dir/Panorama.cpp.i

CMakeFiles/Panorama.dir/Panorama.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Panorama.dir/Panorama.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial/Panorama.cpp -o CMakeFiles/Panorama.dir/Panorama.cpp.s

# Object files for target Panorama
Panorama_OBJECTS = \
"CMakeFiles/Panorama.dir/Panorama.cpp.o"

# External object files for target Panorama
Panorama_EXTERNAL_OBJECTS =

Panorama: CMakeFiles/Panorama.dir/Panorama.cpp.o
Panorama: CMakeFiles/Panorama.dir/build.make
Panorama: /usr/lib/libQt5OpenGL.so.5.15.1
Panorama: /usr/lib/libGL.so
Panorama: /usr/lib/libGLU.so
Panorama: /usr/lib/libQt5Widgets.so.5.15.1
Panorama: /usr/lib/libQt5Gui.so.5.15.1
Panorama: /usr/lib/libQt5Core.so.5.15.1
Panorama: CMakeFiles/Panorama.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Panorama"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Panorama.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Panorama.dir/build: Panorama

.PHONY : CMakeFiles/Panorama.dir/build

CMakeFiles/Panorama.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Panorama.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Panorama.dir/clean

CMakeFiles/Panorama.dir/depend:
	cd /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial /home/ilyes/Documents/3A/MVA/3D_Computer_vision/tp1/Panorama_Initial/CMakeFiles/Panorama.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Panorama.dir/depend

