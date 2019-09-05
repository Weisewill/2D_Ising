ising : 2Dising.o Ising_func.o
	g++ -o ising 2Dising.o Ising_func.o

2Dising.o : 2Dising.cpp Ising_func.h
	g++ -c 2Dising.cpp

Ising_func.o : Ising_func.cpp Ising_func.h
	g++ -c Ising_func.cpp

clean : 
	rm ising 2Dising.o Ising_func.o

