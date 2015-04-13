/* openMP of Gauss-Seidel smoother -by Fang Fang */
#include <omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
/* timing in util.h requires -lrt flag to compile */
#define ff() printf("here is fine!\n");

double norm_res(double *redu, double *blacku, int n, double *A); //calculate norm of residual |Au-f|
void GS(double *redu, double *blacku, double *newredu, double *newblacku, int n, double *A);

static int n_iter = 0;

int main(int argc, char **argv){

    /***** initialize *****/

    if (argc != 2) {
        fprintf(stderr, "Functions need one input as number of discretization!\n");
        abort();
    }

    double *redu, *blacku, *newredu, *newblacku, tic, toc, ini_res, res;

    int n;
    n = atoi(argv[1]);
    double h = 1. / n;

    /* A[0] is diagonal of matrix A 
        A[1] is off-diagonal of matrix A */
    double A[] = {2. / pow(h,2),    -1. / pow(h,2)}; 

    redu = (double *) calloc(n-n/2, sizeof(double)); // for even nodes
    blacku = (double *) calloc(n/2, sizeof(double)); // for odd nodes
    newredu = (double *) calloc(n-n/2, sizeof(double)); // for even nodes
    newblacku = (double *) calloc(n/2, sizeof(double)); // for odd nodes

    ini_res = norm_res(redu, blacku, n, A);  //initial residual
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

    /***** do GS smoother *****/
    tic = omp_get_wtime();
    while ((res > ini_res * 1e-6) && (n_iter < iter_term)){ // terminate after iter_term iterations or res/ini_res < 1e-6
        GS(redu, blacku, newredu, newblacku, n, A);
        if (n_iter % step == 0){ 
            res = norm_res(redu, blacku, n, A);
            printf("After %d iterations, the norm of residule is %e, decreased %e \n", n_iter, res, res/ini_res);
        }
    }
    toc = omp_get_wtime();

    /***** end program *****/
    printf("Time elapsed is %f seconds.\n", toc - tic);
    free(redu);
    free(newredu);
    free(blacku);
    free(newblacku);
    return 0;
}

double norm_res(double *redu, double *blacku, int n, double *A){
    double *res;
    double *u;
    int i;
    u = (double *) malloc(sizeof(double) * n);
    res = (double *) malloc(sizeof(double) * n);

#pragma omp parallel for default(none) shared(redu, u, n) 
    for (i=0; i<n-n/2; ++i) u[2*i] = redu[i];

#pragma omp parallel for default(none) shared(blacku, u, n) 
    for (i=0; i<n/2; ++i) u[2*i+1] = blacku[i];

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
    free(u);
    return n_res;
}

void GS(double *redu, double *blacku, double *newredu, double *newblacku, int n, double *A){
    int i;
    ++n_iter;

    /* if n is even */
    if (n%2 == 0){
/* parallel calculation of newu */
#pragma omp parallel for default(none) shared(n, blacku, newredu, A) 
    for (i = 1; i < n/2; ++i) {
        newredu[i] = (1 - (A[1]*blacku[i-1]+ A[1]*blacku[i]))/A[0];
    }
    newredu[0] = (1 -  A[1]*blacku[0]) / A[0];
#pragma omp barrier

#pragma omp parallel for default(none) shared(n, redu, newblacku, A) 
    for (i = 0; i < n/2-1; ++i) {
        newblacku[i] = (1 - (A[1]*redu[i]+ A[1]*redu[i+1]))/A[0];
    }
    newblacku[n/2] = (1 - A[1]*redu[i])/A[0];
#pragma omp barrier

} /* end of if (n%2 == 0) */

    /* if n is odd */
    if (n%2 == 1){
/* parallel calculation of newu */
#pragma omp parallel for default(none) shared(n, blacku, newredu, A) 
    for (i = 1; i < n/2; ++i) {
        newredu[i] = (1 - (A[1]*blacku[i-1]+ A[1]*blacku[i]))/A[0];
    }
    newredu[0] = (1 -  A[1]*blacku[0]) / A[0];
    newredu[n-n/2] = (1 -  A[1]*blacku[n/2]) / A[0];
#pragma omp barrier

#pragma omp parallel for default(none) shared(n, redu, newblacku, A) 
    for (i = 0; i < n/2-1; ++i) {
        newblacku[i] = (1 - (A[1]*redu[i]+ A[1]*redu[i+1]))/A[0];
    }
#pragma omp barrier

} /* end of if (n%2 == 1) */
    
/* parallel updates of u */
#pragma omp parallel for default(none) shared(n, redu, newredu) 
    for (i = 0; i < n-n/2; ++i) redu[i] = newredu[i];
#pragma omp parallel for default(none) shared(n, blacku, newblacku) 
    for (i = 0; i < n/2; ++i) blacku[i] = newblacku[i];
}

