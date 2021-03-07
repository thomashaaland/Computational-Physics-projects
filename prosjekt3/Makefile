PROG = gauleg
CXX = g++
OBJS = gauleg.o main.o
CXXFLAGS = -fopenmp -Wall 
OUT = Gauss_Legendre_table.txt

gauleg : $(OBJS)
	$(CXX) -o  $(PROG) $(CXXFLAGS) $(OBJS)
main.o :
	$(CXX) -c main.cpp
gauleg.o :
	$(CXX) -c gauleg.cpp
clean:
	rm -f core $(PROG) $(OBJS)
