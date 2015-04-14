CC=g++-4.9
CFLAGS=-std=c++11
INCLUDEBOOST=-I /home/lao/boost_1_57_0
INCLUDEBCDPP=-I /home/lao/Documents/cdboost/include
MODELSHEADERS=atomic-models/reaction.hpp atomic-models/filter.hpp atomic-models/space.hpp
TINYHEADERS=tinyXML/tinyxml.h tinyXML/tinystr.h
STRUCTUREHEADERS=data-structures/types.hpp data-structures/randomNumbers.hpp data-structures/unit_definition.hpp
PARSERHEADERS=parser/parser.hpp

all: main-ammp.o tinyXML/tinyxml.o tinyXML/tinyxmlerror.o tinyXML/tinyxmlparser.o tinyXML/tinystr.o data-structures/unit_definition.o parser/parser.o
	$(CC) -o ammp main-ammp.o tinyXML/tinyxml.o tinyXML/tinyxmlerror.o tinyXML/tinyxmlparser.o tinyXML/tinystr.o data-structures/unit_definition.o parser/parser.o

main-ammp.o: main-ammp.cpp $(MODELSHEADERS) $(TINYHEADERS) $(STRUCTUREHEADERS)
	$(CC) -c $(CFLAGS) $(INCLUDEBOOST) $(INCLUDEBCDPP) main-ammp.cpp -o main-ammp.o

tinyXML/tinyxml.o: tinyXML/tinyxml.cpp $(TINYHEADERS)
	$(CC) -c tinyXML/tinyxml.cpp -o tinyXML/tinyxml.o

tinyXML/tinyxmlerror.o: tinyXML/tinyxmlerror.cpp $(TINYHEADERS)
	$(CC) -c tinyXML/tinyxmlerror.cpp -o tinyXML/tinyxmlerror.o

tinyXML/tinyxmlparser.o: tinyXML/tinyxmlparser.cpp $(TINYHEADERS)
	$(CC) -c tinyXML/tinyxmlparser.cpp -o tinyXML/tinyxmlparser.o

tinyXML/tinystr.o: tinyXML/tinystr.cpp tinyXML/tinystr.h $(TINYHEADERS)
	$(CC) -c tinyXML/tinystr.cpp -o tinyXML/tinystr.o

data-structures/unit_definition.o: data-structures/unit_definition.cpp data-structures/unit_definition.hpp
	$(CC) -c data-structures/unit_definition.cpp -o data-structures/unit_definition.o

parser/parser.o: parser/parser.cpp parser/parser.hpp $(STRUCTUREHEADERS)
	$(CC) -c $(CFLAGS) parser/parser.cpp -o parser/parser.o

clean:
	rm -f ammp *.o *~
	-for d in tinyXML; do (cd $$d; rm -f *.o *~ ); done
	-for d in data-structures; do (cd $$d; rm -f *.o *~ ); done
	-for d in atomic-models; do (cd $$d; rm -f *.o *~ ); done
	-for d in parser; do (cd $$d; rm -f *.o *~ ); done