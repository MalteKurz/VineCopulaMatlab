
# Get the current directory and setting a compiler option for echo
ifeq ($(OS),Windows_NT)
	currentdir = pwd -W
	ECHO = echo -e
	SharedLibName = libVineCopulaCPP.dll
	OperatingSystem = Windows
else
	currentdir = pwd
	ECHO = echo
	SharedLibName = libVineCopulaCPP.so.1.0
	OperatingSystem = Linux
endif

# Directories for Matlab compilation
mingwroot = "C:\\MinGW\\64"
mingw_mex_xml = "mex_C++_mingw-w64.xml"

# Directories
export prefix = /usr
export exec_prefix = $(prefix)

# Include and library directories
export includedir = $(prefix)/include
export libdir = $(exec_prefix)/lib
export oldincludedir = $(includedir)

# Get directories for dependend software
# Boost libraries
export boostdir = $(libdir)
export boostincludedir = $(boostdir)/include
export boostlibdir = $(boostdir)/lib

# NLopt
export nloptdir = $(prefix)
export nloptincludedir = $(nloptdir)/include
export nloptlibdir = $(nloptdir)/lib

INCLUDE_DIRS := -I$(boostincludedir) -I$(nloptlibdir) -I$(includedir)
LIB_DIRS := -L$(nloptlibdir) -L$(libdir)/


CC = g++
CFLAGS = -fPIC -g -c -Wall
FOPEN = -fopenmp

######################################################
LIBRARY_DIR := src
MATLAB_PRIVATE_DIR := private

# Get the list of cpp files
CPP_SRC_FILES := $(wildcard $(MATLAB_PRIVATE_DIR)/*.cpp)
CPP_SRC_FILES := $(filter-out $(MATLAB_PRIVATE_DIR)/SDtau.cpp, $(CPP_SRC_FILES))
DOT_MEX_FILES := $(patsubst %.cpp,%.mexa64,$(CPP_SRC_FILES))

all: $(LIBRARY_DIR)/Makefile CompileLibrary

$(LIBRARY_DIR)/Makefile:
	if [ -d ".git" ]; then	\
	git submodule init;	\
	git submodule update;	\
	else	\
	wget -O VineCopulaCPP.zip https://github.com/MalteKurz/VineCopulaCPP/archive/master.zip;	\
	unzip VineCopulaCPP.zip -d $(LIBRARY_DIR);	\
	mv $(LIBRARY_DIR)/VineCopulaCPP-master/* $(LIBRARY_DIR)/;	\
	rm -r $(LIBRARY_DIR)/VineCopulaCPP-master/;	\
	rm VineCopulaCPP.zip;	\
	fi

CompileLibrary:
	$(MAKE) -C $(LIBRARY_DIR)

InstallLibrary:
	$(MAKE) install -C $(LIBRARY_DIR)

install: InstallLibrary $(OperatingSystem)CompileMatlabMex

WindowsCompileMatlabMex:
	matlab -wait -nodisplay -nodesktop -nosplash -r "mingwroot = '`echo $(mingwroot)`';mingw_mex_xml = '`echo $(mingw_mex_xml)`';libdir = '`echo $(libdir)`';IncludeDirs = '`echo ${INCLUDE_DIRS}`';LibDirs = '`echo ${LIB_DIRS}`';cd('`$(currentdir)`/private'); installVineCopulaMatlab(libdir,IncludeDirs,LibDirs,mingwroot,mingw_mex_xml); exit"

LinuxCompileMatlabMex:
	matlab -nodesktop -nosplash -r "libdir = '`echo $(libdir)`';IncludeDirs = '`echo ${INCLUDE_DIRS}`';LibDirs = '`echo ${LIB_DIRS}`';cd('`$(currentdir)`/private'); installVineCopulaMatlab(libdir,IncludeDirs,LibDirs); exit"

