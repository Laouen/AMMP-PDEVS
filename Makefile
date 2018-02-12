CC=time g++
CFLAGS=-std=c++17
INCLUDE_CADMIUM=-I vendor/cadmium/include
INCLUDE_DESTIME=-I vendor/DESTimes/include
INCLUDE_EXPORTER=-I vendor/DEVSDiagrammer/model_json_exporter/include
INCLUDE_MEMORE=-I vendor/MeMoRe/include
INCLUDE_TINY=-I vendor/tinyxml2
INCLUDE_MONGOCXX = $(shell pkg-config --cflags --libs libmongocxx)
INCLUDE_VENDORS=$(INCLUDE_CADMIUM) $(INCLUDE_DESTIME) $(INCLUDE_EXPORTER) $(INCLUDE_TINY) $(INCLUDE_MEMORE)

INCLUDE_PMGBP=-I include



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

all: check_dirs build/main.o build/types.o build/space.o build/recorder.o build/sink.o vendor/tinyxml2/tinyxml2.o
	$(CC) -g $(CFLAGS) $(INCLUDE_VENDORS) build/main.o build/types.o build/space.o build/recorder.o build/sink.o vendor/tinyxml2/tinyxml2.o -o bin/model $(INCLUDE_MONGOCXX)

build/main.o: check_dirs main.cpp
	$(CC) -g -c $(D) $(CFLAGS) $(INCLUDE_VENDORS) $(INCLUDE_PMGBP) main.cpp -o build/main.o $(shell pkg-config --cflags --libs libmongocxx)

build/space.o: check_dirs src/pmgbp/structures/space.cpp include/pmgbp/structures/space.hpp
	$(CC) -g -c $(D) $(CFLAGS) $(INCLUDE_CADMIUM) $(INCLUDE_PMGBP) src/pmgbp/structures/space.cpp -o build/space.o

build/types.o: check_dirs src/pmgbp/structures/types.cpp include/pmgbp/structures/types.hpp
	$(CC) -g -c $(D) $(CFLAGS) $(INCLUDE_CADMIUM) $(INCLUDE_PMGBP) src/pmgbp/structures/types.cpp -o build/types.o

build/recorder.o: check_dirs vendor/MeMoRe/src/recorder.cpp
	$(CC) -g -c $(CFLAGS) $(INCLUDE_MEMORE) vendor/MeMoRe/src/recorder.cpp -o build/recorder.o $(INCLUDE_MONGOCXX)

build/sink.o: check_dirs vendor/MeMoRe/src/sink.cpp
	$(CC) -g -c $(CFLAGS) $(INCLUDE_MEMORE) vendor/MeMoRe/src/sink.cpp -o build/sink.o $(INCLUDE_MONGOCXX)

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
	rm -rf bin/* build/* *.o *~
