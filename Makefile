EXECUTABLES = omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6 jacobi-omp gs-omp
COMPILER = gcc 
FLAGS = -fopenmp -lm #need to link math library here

all: $(EXECUTABLES)

omp_solved2: omp_solved2.c
	$(COMPILER) $(FLAGS) omp_solved2.c -o omp_solved2 

omp_solved3: omp_solved3.c
	$(COMPILER) $(FLAGS) omp_solved3.c -o omp_solved3 

omp_solved4: omp_solved4.c
	$(COMPILER) $(FLAGS) omp_solved4.c -o omp_solved4 

omp_solved5: omp_solved5.c
	$(COMPILER) $(FLAGS) omp_solved5.c -o omp_solved5 

omp_solved6: omp_solved6.c
	$(COMPILER) $(FLAGS) omp_solved6.c -o omp_solved6 

jacobi-omp: jacobi-omp.c
	$(COMPILER) $(FLAGS) jacobi-omp.c -o jacobi-omp 

gs-omp: gs-omp.c
	$(COMPILER) $(FLAGS) gs-omp.c -o gs-omp 

clean:
	rm -rf $(EXECUTABLES)
