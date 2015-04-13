/* openMP of Jacobi smoother -by Fang Fang */
#include <omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/* timing in util.h requires -lrt flag to compile */
#define ff() printf("here is fine!\n");

double norm_res(double *u, int n, double *A); //calculate norm of residual |Au-f|
void Jacobi(double *u, double *newu, int n, double *A);

static int n_iter = 0;

int main(int argc, char **argv){

    /***** initialize *****/

    if (argc != 2) {
        fprintf(stderr, "Functions need one input as number of discretization!\n");
        abort();
    }

    double *u, *newu, tic, toc, ini_res, res;

    int n;
    n = atoi(argv[1]);
    double h = 1. / n;

    /* A[0] is diagonal of matrix A 
        A[1] is off-diagonal of matrix A */
    double A[] = {2. / pow(h,2),    -1. / pow(h,2)}; 

    u = (double *) calloc(n, sizeof(double)); 
    newu = (double *) calloc(n, sizeof(double)); 

    ini_res = norm_res(u, n, A);  //initial residual
    res = ini_res;

    int iter_term = 1e5;
    static int step = 1e3;
    iter_term = 1e10/n; step = iter_term/10;;
    printf("Program will terminate after %d iterations.\n", iter_term);
    printf("Residual presented every %d iterations.\n", step);
  
    #pragma omp parallel 
    {
     if (omp_get_thread_num() == 0) {
         printf("Number of threads = %d\n", omp_get_num_threads()); }
    } 

    /***** do Jacobi smoother *****/
    tic = omp_get_wtime();
    while ((res > ini_res * 1e-6) && (n_iter < iter_term)){ // terminate after iter_term iterations or res/ini_res < 1e-6
        Jacobi(u, newu, n, A);
        if (n_iter % step == 0){ 
            res = norm_res(u, n, A);
            printf("After %d iterations, the norm of residule is %e, decreased %e \n", n_iter, res, res/ini_res);
        }
    }
    toc = omp_get_wtime();

    /***** end program *****/
    printf("Time elapsed is %f seconds.\n", toc - tic);
    free(u);
    free(newu);
    return 0;
}

double norm_res(double *u, int n, double *A){
    double *res;
    int i;
    res = (double *) malloc(sizeof(double) * n);

/* parallel calculation of residual */
    res[0] = A[0]*u[0] + A[1]*u[1]  - 1;
#pragma omp parallel for default(none) shared(res, u, n, A) 
    for (i = 1; i < n-1; i ++){
    res[i] =  A[1]*u[i-1] + A[0]*u[i] + A[1]*u[i+1] - 1;
    }
    res[n-1] = A[1]*u[n-2] + A[0]*u[n-1] - 1;
    
/* parallel sum up */
    double n_res = 0.;
#pragma omp parallel for reduction(+:n_res) 
    for (i = 0; i < n; ++i){
        n_res += res[i]*res[i];
    }
    n_res = sqrt(n_res);

    free(res);
    return n_res;
}

void Jacobi(double *u, double *newu, int n, double *A){
    int i;
    ++n_iter;

/* parallel calculation of newu */
    newu[0] = (1 -  A[1]*u[1]) / A[0];
#pragma omp parallel for default(none) shared(n, u, newu, A) 
    for (i = 1; i < n-1; ++i) {
        newu[i] = (1 - (A[1]*u[i-1]+ A[1]*u[i+1]))/A[0];
    }
    newu[n-1] = (1 - A[1]*u[n-2])/A[0];

/* parallel updates of u */
#pragma omp parallel for default(none) shared(n, u, newu) 
    for (i = 0; i < n; ++i) {
        u[i] = newu[i];
    }
}

