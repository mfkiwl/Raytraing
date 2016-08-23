main: msis00.o main.o nrlmsise00_sub.o
	gfortran -o main main.o msis00.o nrlmsise00_sub.o 
main.o: main.c
	cc -Wall -c main.c
msis00.o: msis00.c msis00.h
	gfortran -Wall -c msis00.c 
nrlmsise00_sub.o: nrlmsise00_sub.for
	gfortran -Wall -c nrlmsise00_sub.for

