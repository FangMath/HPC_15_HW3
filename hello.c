/*
 * OpenMP Example based on code in Tim Mattson's Tutorial
 *
 * [slides](http://openmp.org/mp-documents/Intro_To_OpenMP_Mattson.pdf)
 * [video](http://goo.gl/EaxWjY)
 */

#include <stdio.h>
#include <omp.h>

int main(int argc, char* argv[])
{
  double tic, toc;
  /* omp_set_num_threads(8); */
  tic = omp_get_wtime();
#pragma omp parallel
  {
    int id = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    printf("hello(%d) ",  id);
    printf("world(%d)\n", id);
    if(id == 0)
      printf("num_threads = %d\n", num_threads);
  }
  toc = omp_get_wtime();

  printf("elapsed time = %g s\n", toc-tic);

  return 0;
}
