
LIBRARY_DIR := src

all:
	git submodule init
	git submodule update
	$(MAKE) -C $(LIBRARY_DIR)
	matlab -nodesktop -nosplash -r "cd('`pwd`/private'); installVineCopulaMatlab; exit"

