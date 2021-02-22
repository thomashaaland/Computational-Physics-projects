OBJS = main.o jacobi.o

jacobi_method : $(OBJS)
	g++ -o jacobi_method $(OBJS)
main.o : 
	g++ -c main.cpp
jacobi.o : jacobi.h
	g++ -c jacobi.cpp
clean:
	rm -f core jacobi_method $(OBJS)
