all: lib/liblocalsearch.a

lib/liblocalsearch.a: build/localsearch.o
	ar rcs lib/liblocalsearch.a build/*.o

build/localsearch.o: src/localsearch.cpp include/liblocalsearch/localsearch.h
	g++ -ansi -I ./include -c src/localsearch.cpp -o build/localsearch.o

clean: .
	rm build/*.o lib/liblocalsearch.a
