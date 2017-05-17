CC=g++
CFLAGS=-std=c++14
INCLUDECADMIUM=-I ../cadmium/include
INCLUDEDESTIME=-I ../DESTimes
STRUCTUREHEADERS=data-structures/types.hpp data-structures/randomNumbers.hpp data-structures/unit_definition.hpp

all: main.o
	$(CC) -g -o model main.o

main.o: main.cpp $(STRUCTUREHEADERS)
	$(CC) -g -c $(CFLAGS) $(INCLUDECADMIUM) $(INCLUDEDESTIME) main.cpp -o main.o

data-structures/unit_definition.o: data-structures/unit_definition.cpp data-structures/unit_definition.hpp
	$(CC) -g -c $(CFLAGS) $(INCLUDEBOOST) $(INCLUDEBCDPP) data-structures/unit_definition.cpp -o data-structures/unit_definition.o

data-structures/types.o: data-structures/types.cpp data-structures/types.hpp
	$(CC) -g -c $(CFLAGS) $(INCLUDEBOOST) $(INCLUDEBCDPP) data-structures/types.cpp -o data-structures/types.o

clean:
	rm -f ammp *.o *~
	-for d in data-structures; do (cd $$d; rm -f *.o *~ ); done
	-for d in atomic-models; do (cd $$d; rm -f *.o *~ ); done