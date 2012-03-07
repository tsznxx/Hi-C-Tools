OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/gzstream/ -I$(UTILITIES_DIR)/stringUtils/ 

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= fileUtils.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=gzstream.o stringUtils.o
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
