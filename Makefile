
LIBRARY_DIR := src
MATLAB_PRIVATE_DIR := private

# Get the list of cpp files
CPP_SRC_FILES := $(wildcard $(MATLAB_PRIVATE_DIR)/*.cpp)
CPP_SRC_FILES := $(filter-out $(MATLAB_PRIVATE_DIR)/SDtau.cpp, $(CPP_SRC_FILES))
DOT_MEX_FILES := $(patsubst %.cpp,%.mexa64,$(CPP_SRC_FILES))

all: $(LIBRARY_DIR)/Makefile CompileLibrary $(DOT_MEX_FILES)

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

$(DOT_MEX_FILES):
	matlab -nodesktop -nosplash -r "cd('`pwd`/private'); installVineCopulaMatlab; exit"

