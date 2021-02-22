#Makefile for jacobi application
#Created by Thomas Haaland 22. february 2021

PROG = jacobi
CXX = g++
CXXFLAGS = -std=c++11 -O2
LDFLAGS = -larmadillo
OBJS = main.o jacobi.o test.o
OUT = eigvecs_u.dat iterations_n.dat solutions.dat

jacobi : $(OBJS)
	$(CXX) -o $(PROG) $(OBJS) $(LDFLAGS) $(CXXFLAGS)
main.o :
	$(CXX) -c $(CXXFLAGS) main.cpp
jacobi.o : jacobi.h
	$(CXX) -c $(CXXFLAGS) jacobi.cpp
test.o : test.h
	$(CXX) -c $(CXXFLAGS) test.cpp
clean:
	rm -f core $(PROG) $(OBJS) $(OUT)
