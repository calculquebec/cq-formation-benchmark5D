EXEC=b5D.omp
OBJ=cubicmatrix.o cmvectorized.o cmopenmp.o
OBJ_OMP=main.o

all: $(EXEC)

b5D.omp: $(OBJ) $(OBJ_OMP)
	g++ -fopenmp -O3 -o $@ $+

%.o: %.cpp %.h cubicmatrix.h
	g++ -c -fopenmp -O3 -o $@ $<

clean:
	rm -f $(EXEC) $(OBJ_OMP) $(OBJ)
