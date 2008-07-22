all: build/liblocalsearch.a

build/liblocalsearch.a: build/localsearch.o
	ar rcs build/liblocalsearch.a build/*.o

build/localsearch.o: src/localsearch.cpp include/liblocalsearch/localsearch.h
	g++ -ansi -I ./include -c src/localsearch.cpp -o build/localsearch.o

doc: ./src/ ./include/liblocalsearch/
	doxygen

clean: .
	rm build/*.o build/liblocalsearch.a
