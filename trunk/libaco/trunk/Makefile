all: build/libaco.a

build/libaco.a: build/ants.o build/graph.o build/util.o
	ar rcs build/libaco.a build/*.o

build/ants.o: src/ants.cpp include/libaco/ants.h include/libaco/graph.h include/libaco/util.h
	g++ -Wall -ansi -I ./include -c src/ants.cpp -o build/ants.o

build/graph.o: src/graph.cpp include/libaco/graph.h
	g++ -Wall -ansi -I ./include -c src/graph.cpp -o build/graph.o

build/util.o: src/util.cpp include/libaco/util.h
	g++ -Wall -ansi -I ./include -c src/util.cpp -o build/util.o

doc: src/ include/libaco/
	doxygen

clean:
	rm build/*.o build/libaco.a
