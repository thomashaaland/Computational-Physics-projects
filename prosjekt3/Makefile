PROG = gauleg
CXX = g++
OBJS = gauleg.o main.o
OUT = Gauss_Legendre_table.txt

gauleg : $(OBJS)
	$(CXX) -o $(PROG) $(OBJS)
main.o :
	$(CXX) -c main.cpp
gauleg.o :
	$(CXX) -c gauleg.cpp
clean:
	rm -f core $(PROG) $(OBJS)
