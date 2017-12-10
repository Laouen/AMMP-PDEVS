CC=g++
CFLAGS=-std=c++17 -ftemplate-depth=900000
INCLUDECADMIUM=-I vendor/cadmium/include
INCLUDEDESTIME=-I vendor/DESTimes/include
INCLUDEEXPORTER=-I vendor/DEVSDiagrammer/model_json_exporter
INCLUDETINY=-I vendor/tinyxml2
VENDOR=$(INCLUDECADMIUM) $(INCLUDEDESTIME) $(INCLUDEEXPORTER) $(INCLUDETINY)
INCLUDE=-I include

# =============== Parameters ==================== #
# D: all the -D flags, they are:
#  * DIAGRAM: If this flag is set, the model will compile the DEVSDiagrammer mode and no simulation will run
#  * show_log: Make the logger to print the log messages. 
#  * show_info: Make the logger to print the info messages. 
#  * show_debug: Make the logger to print the debug messages. 
#  * show_error: Make the logger to print the error messages. 
#  
#  example: D='-D DIAGRAM' will compile the model in the DEVSDiagrammer mode and the model diagram .json will be print
# ================================================ #

all: build/main.o
	$(CC) -g build/main.o build/types.o build/space.o vendor/tinyxml2/tinyxml2.o -o bin/model

build/main.o: main.cpp build/types.o build/space.o vendor/tinyxml2/tinyxml2.o
	$(CC) -g -c $(D) $(CFLAGS) $(VENDOR) $(INCLUDE) main.cpp -o build/main.o

build/space.o: check_dirs src/structures/space.cpp include/pmgbp/structures/space.hpp
	$(CC) -g -c $(D) $(CFLAGS) $(INCLUDECADMIUM) $(INCLUDE) src/structures/space.cpp -o build/space.o

build/types.o: check_dirs src/structures/types.cpp include/pmgbp/structures/types.hpp
	$(CC) -g -c $(D) $(CFLAGS) $(INCLUDECADMIUM) $(INCLUDE) src/structures/types.cpp -o build/types.o

vendor/tinyxml2/tinyxml2.o: vendor/tinyxml2/tinyxml2.h vendor/tinyxml2/tinyxml2.cpp
	$(CC) -g -c $(CFLAGS) vendor/tinyxml2/tinyxml2.cpp -o vendor/tinyxml2/tinyxml2.o

.PHONY: clean clean_model clean_all check_dirs

check_dirs:
	mkdir -p bin
	mkdir -p build

clean_all: clean_model clean

clean_model:
	rm top.hpp top_ports.hpp parameters.xml

clean:
	rm -rf bin/ build/ *.o *~
