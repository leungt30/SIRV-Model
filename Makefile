all: sir

sir: sir.o
	c++ -o sir sir.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

sir.o: sir.cpp
	c++ -c sir.cpp
