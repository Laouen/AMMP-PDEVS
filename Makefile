CC=time g++
CFLAGS=-std=c++14 -ftemplate-depth=900000
INCLUDECADMIUM=-I ../cadmium/include
INCLUDEDESTIME=-I ../DESTimes/include
INCLUDEEXPORTER=-I vendor/DEVSDiagrammer/model_json_exporter
INCLUDETINY=-I vendor/tinyxml2
INCLUDELIBS=-I libs
INCLUDESTRUCTURES=-I structures

all: main.o
	$(CC) -g -o model main.o structures/types.o structures/space.o vendor/tinyxml2/tinyxml2.o

main.o: main.cpp structures/types.o structures/space.o vendor/tinyxml2/tinyxml2.o
	$(CC) -g -c $(CFLAGS) $(INCLUDECADMIUM) $(INCLUDEDESTIME) $(INCLUDEEXPORTER) $(INCLUDETINY) $(INCLUDELIBS) $(INCLUDESTRUCTURES) main.cpp -o main.o

structures/space.o: structures/space.cpp structures/space.hpp
	$(CC) -g -c $(CFLAGS) $(INCLUDECADMIUM) $(INCLUDELIBS) structures/space.cpp -o structures/space.o

structures/types.o: structures/types.cpp structures/types.hpp
	$(CC) -g -c $(CFLAGS) $(INCLUDECADMIUM) $(INCLUDELIBS) structures/types.cpp -o structures/types.o

vendor/tinyxml2/tinyxml2.o: vendor/tinyxml2/tinyxml2.h vendor/tinyxml2/tinyxml2.cpp
	$(CC) -g -c $(CFLAGS) vendor/tinyxml2/tinyxml2.cpp -o vendor/tinyxml2/tinyxml2.o

clean:
	rm -f model *.o *~
	-for d in structures; do (cd $$d; rm -f *.o *~ ); done
	-for d in vendor/tinyxml2; do (cd $$d; rm -f *.o *~ ); done