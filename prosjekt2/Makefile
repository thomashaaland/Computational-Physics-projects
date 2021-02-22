#Makefile for jacobi application
#Created by Thomas Haaland 22. february 2021

PROG = jacobi
CC = g++
CPPFLAGS = -std=c++11 -O2
LDFLAGS = -larmadillo
OBJS = main.o jacobi.o test.o
OUT = eigvecs_u.dat iterations_n.dat solutions.dat

jacobi : $(OBJS)
	$(CC) -o $(PROG) $(LDFLAGS) $(OBJS)
main.o :
	$(CC) -c $(CPPFLAGS) main.cpp
jacobi.o : jacobi.h
	$(CC) -c $(CPPFLAGS) jacobi.cpp
test.o : test.h
	$(CC) -c $(CPPFLAGS) test.cpp
clean:
	rm -f core $(PROG) $(OBJS) $(OUT)
