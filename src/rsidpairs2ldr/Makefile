UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

CXXFLAGS = -Wall -O3 -D_FILE_OFFSET_BITS=64 -fPIC -lz

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/common/ \
		   -I$(UTILITIES_DIR)/stringUtils/ \
           -I$(UTILITIES_DIR)/gzstream/ \

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= rsidpairs2ldr.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=gzstream.o 
PROGRAM= rsidpairs2ldr

EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
BUILT_PROGRAM= $(patsubst %,$(BIN_DIR)/%,$(PROGRAM))



all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@echo $(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)
	@echo $(CXX) -o $(BUILT_PROGRAM) $(BUILT_OBJECTS) $(EXT_OBJECTS) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)
	@$(CXX) -o $(BUILT_PROGRAM) $(BUILT_OBJECTS) $(EXT_OBJECTS) $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(ALL_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/lineFileUtilities/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bedUtils/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/genomeFile/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/gzstream/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/fileType/	
		
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
