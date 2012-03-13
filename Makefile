# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

SHELL := /bin/bash -e

# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export UTIL_DIR	= src/utils
export CXX		= g++
export CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC
export LIBS		= -lz
export BT_ROOT  = src/utils/BamTools/


SUBDIRS = $(SRC_DIR)/test \
	      $(SRC_DIR)/genomeToREFrags \
		  $(SRC_DIR)/endsMappAbility \
#		  $(SRC_DIR)/bedsAnnotation
#		  $(SRC_DIR)/intersectBed 

UTIL_SUBDIRS =	$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/BamTools \
				$(SRC_DIR)/utils/seqUtils \
				$(SRC_DIR)/utils/fileUtils \
				$(SRC_DIR)/utils/bedUtils

BUILT_OBJECTS = $(OBJ_DIR)/*.o

all:
	[ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)
	
	@echo "Building BEDTools:"
	@echo "========================================================="
	
	@for dir in $(UTIL_SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

#	@echo "- Building main bedtools binary."
#	@$(CXX) $(CXXFLAGS) -c src/bedtools.cpp -o obj/bedtools.o -I$(UTIL_DIR)/version/
#	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/bedtools $(BUILT_OBJECTS) -L$(UTIL_DIR)/BamTools/lib/ -lbamtools $(LIBS)
#	@echo "done."
#	
#	@echo "- Creating executables for old CLI."
#	@python scripts/makeBashScripts.py
#	@chmod +x bin/*
#	@echo "done."
	

.PHONY: all

clean:
	@echo "Cleaning up."	
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
	@rm -Rf $(BT_ROOT)/lib
	@rm -f $(BT_ROOT)/src/api/*.o
	@rm -f $(BT_ROOT)/src/api/internal/*.o
	@rm -Rf $(BT_ROOT)/include

install:
	@cp bin/* $(HOME)/bin

.PHONY: clean
