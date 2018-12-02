# -------------------------------------------------------------------------
# To compile tools, cd to directory with this 'makefile' and use
#
# make ExtractMPM
# make CompareGlobal
#
# to compile a tool. Use make along to compile all tools.
# Each make command can be modified with the following options
#
# Add 'SYSTEM=system' where the lower case system is where
#    you are compiling. Supported systems are:
#       carbon - Unix cluster at OSU
#       mac - MacOS X using gcc
#       mac-clang - MacOS X using clang-mp
#
# Other options are
#    make clean - to remove all compiled tool objects and executables
#    make install - to copy all compiled tools to desired installation folder
#
# If you need other options, you can edit the makefile or override them in the make command
# The following are the most important variables
#
# 1. $(CC) is the gcc compiler, $(RM) is the rm command, $(CP) is copy for install
# 2. $(CFLAGS) is flags for gcc compiler options
#    $(LFLAGS) is flags fpr gcc linking options
# 3. $(ioutput) is path to install folder (default is "~/bin")
#
# Each of these can changed by editing or can be overridden at make time as
# documented more in the numbered section below
#
# -------------------------------------------------------------------------

# defatul to build all
.PHONY: default
default: all

# Default system - OSU Linux cluster (old name, but applied to current)
SYSTEM=carbon

# 1. Define compiler to use and other commands
#     CC is path to the c++ compiler, must support OpenMP is used
#     RM is path/options to the file remove command (used only in make clean)
#     CP is path/options to copy command (used only in make install
CC= g++
ifeq ($(SYSTEM),carbon)
    CC= /usr/local/gcc48/bin/g++
endif
ifeq ($(SYSTEM),mac)
    CC= /opt/local/bin/g++-mp-4.9
endif
ifeq ($(SYSTEM),mac-clang)
    CC= /opt/local/bin/clang-mp-4.0
endif
RM= rm
CP= cp

# 2. Compiler options -c=compile only is required, others are optional
#     CFLAGS are flags passed to compile each source file. If must contain -c and
#	  can add any other needed options. Some flags found useful in this project are
#            -O3 = optimize most -O0 = not optimized
#            -fopenmp = enable OpenMP compile (link too)
#            -g = debugging info  -pg = profiling
#            -arch i386 = compile architecture (MacOS X only, other archs are x86_64, ppc, and ppc64)
#     LFLAGS are flags passed at link time. The default is -fopenmp. Some flags found useful
#       in this project are
#            -arch i386 = link architecture (MacOS X only, other archs are x86_64, ppc, and ppc64)
#            -pg = profiling
#     If using -arch, or -pg, must have in both CFLAGS and LFLAGS. For MacOS X, the specified arch
#       must match the xerces library used for linking
CFLAGS= -c -O3
LFLAGS=
ifeq ($(SYSTEM),mac-clang)
    LFLAGS= -lc++
endif

# 3. Define paths to intall folder (relative to 'makefile')
ioutput = ~/bin

# -------------------------------------------------------------------------
# all compiled objects and executables
objects = ExtractMPM.o CompareGlobal.o
tools = ExtractMPM CompareGlobal

# -------------------------------------------------------------------------
# Default to make ExtractMPM tool
.PHONY : all
all : ExtractMPM CompareGlobal

# -------------------------------------------------------------------------
# Tool targets
ExtractMPM : ExtractMPM.o
	$(CC) $(LFLAGS) -o ExtractMPM ExtractMPM.o
CompareGlobal : CompareGlobal.o
	$(CC) $(LFLAGS) -o CompareGlobal CompareGlobal.o

# -------------------------------------------------------------------------
# Compile
ExtractMPM.o : ExtractMPM.cpp ExtractMPM.hpp
	$(CC) $(CFLAGS) -o ExtractMPM.o ExtractMPM.cpp
CompareGlobal.o : CompareGlobal.cpp CompareGlobal.hpp
	$(CC) $(CFLAGS) -o CompareGlobal.o CompareGlobal.cpp

# -------------------------------------------------------------------------
# To clean compiled objects        
clean : 
	$(RM) $(objects) $(tools)

# -------------------------------------------------------------------------
# To install        
install : 
	$(CP) $(tools) $(ioutput)
