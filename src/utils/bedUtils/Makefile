OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/common/ -I$(UTILITIES_DIR)/gzstream/ -I$(UTILITIES_DIR)/stringUtils/ -I$(UTILITIES_DIR)/seqUtils/ -I$(UTILITIES_DIR)/fileUtils/ 

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= bedUtils.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=gzstream.o seqUtils.o fileUtils.o 
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C -W $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
