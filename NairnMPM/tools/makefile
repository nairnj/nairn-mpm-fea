# -------------------------------------------------------------------------
# To compile ExtractMPM use the following commands
#
# cd to directory with this 'makefile'
# give the command 'make SYSTEM=system' where the lower case system is where
#    you are compiling. Supported systems are:
#       mac - MacOS X using gcc
#       mac-clang - MacOS X using clang-mp
#       mpm - OSU machine
#       cygwin - Standard cygwin install in Windows
#		cygwinbin - Standard cygwin install using binary xerces
# Other options are
#    make clean - to remove all compiled objects and executable
#    make mostlyclean - to remove all compiled objects only
#    make install - to copy compiled code to desired installation destination
#
# If you need other options, you can edit the makefile or override them in the make command
# The following are the most important variables
#
# 1. $(CC) is the gcc compiler, $(RM) is the rm command, $(CP) is copy for install
# 2. $(CFLAGS) is flags for gcc compiler options
#    $(LFLAGS) is flags fpr gcc linking options
# 3. $(output) is path to save executable
# 4. $(ioutput) is path to install executable
#
# Each of these can changed by editing or can be overridden at make time as
# documented more in the numbered section below
#
# -------------------------------------------------------------------------
#
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

# 3. Define executable destination
#     output = relative or full path to save executable or to find it for a subsequent make install
#     dtdpath = relative of full path to DTD file, only used during make install
output = ExtractMPM

# 4. Define paths to all files (relative to 'makefile')
#     These are usually not changed unless you reorganize the code
ioutput = ~/bin

# -------------------------------------------------------------------------
# all compiled objects
objects = ExtractMPM.o

# -------------------------------------------------------------------------
# Link all objects
ExtractMPM : $(objects)
	$(CC) $(LFLAGS) -o $(output) $(objects)

# -------------------------------------------------------------------------
# Compile
ExtractMPM.o : ExtractMPM.cpp ExtractMPM.hpp
	$(CC) $(CFLAGS) -o ExtractMPM.o ExtractMPM.cpp

# -------------------------------------------------------------------------
# To make executable
all : NairnMPM
        
# -------------------------------------------------------------------------
# To clean compiled objects        
clean : 
	$(RM) $(objects) $(output)

# -------------------------------------------------------------------------
# To clean compiled objects only      
mostlyclean :
	$(RM) $(objects)

# -------------------------------------------------------------------------
# To install        
install : 
	$(CP) $(output) $(ioutput)
