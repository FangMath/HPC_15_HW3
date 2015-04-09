EXECUTABLES = hello omp_solved2 omp_solved3 #omp_solved4 omp_solved5 omp_solved6
COMPILER = gcc 
FLAGS = -fopenmp

all: $(EXECUTABLES)

hello: hello.c
	$(COMPILER) $(FLAGS) hello.c -o hello 

omp_solved2: omp_solved2.c
	$(COMPILER) $(FLAGS) omp_solved2.c -o omp_solved2 

omp_solved3: omp_solved3.c
	$(COMPILER) $(FLAGS) omp_solved3.c -o omp_solved3 

#omp_solved4: omp_solved4.c
#	$(COMPILER) $(FLAGS) omp_solved4.c -o omp_solved4 
#
#omp_solved5: omp_solved5.c
#	$(COMPILER) $(FLAGS) omp_solved5.c -o omp_solved5 
#
#omp_solved6: omp_solved6.c
#	$(COMPILER) $(FLAGS) omp_solved6.c -o omp_solved6 
#
clean:
	rm -rf $(EXECUTABLES)
