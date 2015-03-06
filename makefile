CC=g++-4.9
CFLAGS=-std=c++11
INCLUDEBOOST=-I /home/lao/boost_1_57_0
INCLUDEBCDPP=-I /home/lao/Documents/cdboost/include
HEADERS=atomic-models/reaction.hpp data-structures/reaction_input.hpp atomic-models/filter.hpp atomic-models/controler.hpp
TINYHEADERS=tinyXML/tinyxml.h tinyXML/tinystr.h

all: main-ammp.o tinyXML/tinyxml.o tinyXML/tinyxmlerror.o tinyXML/tinyxmlparser.o tinyXML/tinystr.o
	$(CC) -o ammp main-ammp.o tinyXML/tinyxml.o tinyXML/tinyxmlerror.o tinyXML/tinyxmlparser.o tinyXML/tinystr.o

main-ammp.o: main-ammp.cpp $(HEADERS) $(TINYHEADERS)
	$(CC) -c $(CFLAGS) $(INCLUDEBOOST) $(INCLUDEBCDPP) main-ammp.cpp -o main-ammp.o

tinyXML/tinyxml.o: tinyXML/tinyxml.cpp $(TINYHEADERS)
	$(CC) -c tinyXML/tinyxml.cpp -o tinyXML/tinyxml.o

tinyXML/tinyxmlerror.o: tinyXML/tinyxmlerror.cpp $(TINYHEADERS)
	$(CC) -c tinyXML/tinyxmlerror.cpp -o tinyXML/tinyxmlerror.o

tinyXML/tinyxmlparser.o: tinyXML/tinyxmlparser.cpp $(TINYHEADERS)
	$(CC) -c tinyXML/tinyxmlparser.cpp -o tinyXML/tinyxmlparser.o

tinyXML/tinystr.o: tinyXML/tinystr.cpp tinyXML/tinystr.h $(TINYHEADERS)
	$(CC) -c tinyXML/tinystr.cpp -o tinyXML/tinystr.o

clean:
	rm -f ammp *.o *~
	-for d in tinyXML; do (cd $$d; rm -f *.o *~ ); done