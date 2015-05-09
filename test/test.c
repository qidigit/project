#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <omp.h>

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif


double f(double x, void* data)
{
	return sin(x);
}

int main(int argc, char* argv[])
{
    int N = 2048;
    int div = 100;
    int howmany = 100000/div;
    int idist = N, odist = N/2+1;
    int istride = 1, ostride = 1;
    int inembed = N, onembed = N; 

    double *u, *u_out;
    fftw_complex *u_hat;
    fftw_plan p1, p2;

    u_hat  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1) * howmany);
    u = (double*) fftw_malloc(sizeof(double) * N * howmany);
    u_out = (double*) fftw_malloc(sizeof(double) * N * howmany);

    fftw_init_threads();
    int nthreads = 8; //omp_get_max_threads();
    fftw_plan_with_nthreads(nthreads);
#pragma omp parallel
    {
      printf("Thread %d of %d initialized...\n", omp_get_thread_num(), omp_get_num_threads());
    }
    

    p1 = fftw_plan_many_dft_r2c(1, &N, howmany,
                               u,     &inembed, istride, idist,
                               u_hat, &onembed, ostride, odist, FFTW_MEASURE);
    p2 = fftw_plan_many_dft_c2r(1, &N, howmany,
                               u_hat, &inembed, istride, odist,
                               u_out, &onembed, ostride, idist, FFTW_MEASURE);
//    p2 = fftw_plan_dft_c2r_1d(N, u_hat, u_out, FFTW_MEASURE);

    int i;
    for (i = 0; i < N*howmany; i++) {
        u[i] = 1 / (1 + (double) i*i/N/N);
    }

    double tic, toc, elapsed;
    tic = omp_get_wtime();
    #pragma omp parallel for private(i)
    for (i=0; i<div; i++) {
    fftw_execute_dft_r2c(p1, u, u_hat);
/*
    for (i = 0; i < N/2+1; i++) {
        printf("u_hat[%d] = %f+%fi,\t%f+%fi,\t%f+%fi\n", i, creal(u_hat[i]), cimag(u_hat[i]),creal(u_hat[i+odist]),cimag(u_hat[i+odist]),creal(u_hat[i+2*odist]), cimag(u_hat[i+2*odist]));
    }

    printf("\n");
    */

    fftw_execute_dft_c2r(p2, u_hat, u_out);
    }
    toc = omp_get_wtime();
    elapsed = toc-tic;
    printf("\nrun_time = %f seconds.\n", elapsed);

    /*
    printf("\n");
    for (i = 0; i < N; i++) {
        printf("u1[%d] = %f,\t%f\n", i, u[i], u_out[i]/N);
    }
    printf("\n");
    for (i = 0; i < N; i++) {
        printf("u2[%d] = %f,\t%f\n", i, u[i+idist], u_out[i+idist]/N);
    }
    printf("\n");
    for (i = 0; i < N; i++) {
        printf("u3[%d] = %f,\t%f\n", i, u[i+2*idist], u_out[i+2*idist]/N);
    }
    */

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(u_hat);
    fftw_free(u);
    fftw_free(u_out);


    return 0;
}

