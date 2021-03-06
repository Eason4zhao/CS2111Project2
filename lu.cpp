#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "lapacke.h"
#include "blas.h"
#include "cblas.h"

using namespace std;

void mydgetrf(double* m, int N, int *pvt, double* tempv) {
        double max;
        int maxind, temp;
        for(int i = 0; i < N-1; i++) {
                maxind = i;
                max = fabs(m[i*N+i]);
                for(int t = i+1; t < N; t++) {
                        if(fabs(m[t*N+i]) > max) {
                                maxind = t;
                                max = fabs(m[t*N+i]);
                        }
                }
        if(max == 0) {
                        printf("LU factoration failed: coefficient matrix is singular.\n");
                        return -1;
                     }
        else {
              if(maxind != i) {
              //save pivoting infomation
              temp = pvt[i];
              pvt[i] = pvt[maxind];
              pvt[maxind] = temp;
              //swap rows
              for(int j = 0; j < N; j++) {
                  // double tempv = m[i*N+j];
                  tempv[j] = m[i*N+j];
                  m[i*N+j] = m[maxind*N+j];
                  m[maxind*N+j] = tempv[j];
              }
            }
              //factorization
              for(int j = i+1; j < N; j++) {
                  m[j*N+i] = m[j*N+i]/m[i*N+i];
                      for(int k = i+1; k < N; k++) {
                          m[j*N+k] = m[j*N+k] - m[j*N+i]*m[i*N+k];
                                }
                        }
                }

        }
        return 0;
}

viod mydtrsm_forward(int N, double* A, int* pvt, double* b, double* x, double *y) {

        double sum = 0.0;
        // double temp;

        // double* y = (double*)calloc(sizeof(double), N);
        y[0] = b[pvt[0]];
        for(int i = 1; i < N; i++) {
                for(int j = 0; j < i-1; j++) {
                        sum += y[j] * A[j*N+j];
                }
                y[i] = b[pvt[i]] - sum;
        }

        return 0;

}


void mydtrsm_back(int N, double* A, int* pvt, double* b, double* x, double* y) {

        double sum = 0.0;
        x[N-1] = y[N-1]/A[N*N-1];
        for(int i = N-2; i >= 0; i--) {
                double sum = 0.0;
                for(int j = i+1; j < N; j++) {
                        sum += x[j]*A[i*N+j];
                }
                x[i] = (y[i] - sum)/A[i*N+i];
        }
        return 0;

}

void transpose(double *a, int N) {
        double temp;
        for(int i=0; i<N; i++) {
                for(int j = i; j < N; j++) {
                        temp = a[i*N+j];
                        a[i*N+j] = a[i*N+i];
                        a[j*N+i] = temp;
                }
        }
}

int main(int argc, char* argv[]) {
        if(argc > 1) {
                int N=0;
                N = atoi(argv[1]);

                int INFO = N;
                int LDA = N;
                int LDB = N;

                int NRHS = 1;
                int *IPIV = (int *)calloc(sizeof(int), N);

                double* B = (double*)calloc(sizeof(double), N);
                double* B1 = (double*)calloc(sizeof(double), N);

                int M = 1;
                double a = 1.0;

                double *A = NULL;
                A  = (double *) malloc(N*N*sizeof(double));

                double *A1 = NULL;
                A1 = (double *) malloc(N*N*sizeof(double));

                double *b = NULL;
                b = (double *) malloc(N*N*sizeof(double));

                double *x = NULL;
                x = (double *) malloc(N*N*sizeof(double));

                double* tempv = NULL;
                tempv = (double*)calloc(sizeof(double), N);

                double max, min;
                max = N;
                min = 0.0;
                double range = 0.0+(max - min);
                double div = RAND_MAX / range;
                srand(time(NULL));

                for(int i = 0; i < N; i++) {
                        B[i] = min + (rand()/div);
                        B1[i] = B[i];
                }

                for(int index = 0; index < N*N; index++) {
                        A[index] = min + (rand()/div);
                        A1[index] = A[index];
                }

                for(int index = 0; index < N; index++) {
                        b[index] = min + (rand()/div);
                }



                int *pvt = (int*)calloc(sizeof(int), N);

                for(int i = 0; i < N; i++) {
                        pvt[i] = i;
                }

                double rtime;
                double gflops;

                double time1 = clock();

                mydgetrf(A, N, pvt, tempv);
                mydtrsm_forward(N, A, pvt, b, x, y)
                mydtrsm_back(N, A, pvt, b, x, y)

                double end1 = clock();

                rtime = (double)end1 - time1/CLOCKS_PER_SEC;
                gflops = (2.0*powf(N,3.0))/(3.0*rtime*powf(10.0,9.0));
                printf("Running time is %.6f seconds.\n", rtime);
                printf("GFLOPS is %f\n", gflops);

                int trans = transpose(A1, N);

                double time2 = clock();

                // LU factorization
                LAPACK_dgetrf(&N,&N,A1,&LDA,IPIV,&INFO);
                cblas_dtrsm(CblasColMajor,CblasLeft,CblasLower,CblasTrans,CblasUnit,N,M,a,A1, N, B, N);
                cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,N,M,a,A1, N, B, N);

                double end2 = clock();

                rtime = (double)end2 - time2/CLOCKS_PER_SEC;
                gflops = (2.0*powf(N,3.0))/(3.0*rtime*powf(10.0,9.0));
                printf("LAPACK Running time is %.6f seconds.\n", rtime);
                printf("LAPACK GFLOPS is: %f\n", gflops);

                trans = transpose(A1, N);

                for(int i = 0; i < N; i++)
                {
                        double tmp = B[IPIV[i]-1];
                        B[IPIV[i]-1] = B[i];
                        B[i] = tmp;
                }

                // forward  L(Ux) = B => y = Ux
                //cblas_dtrsm(CblasColMajor,CblasLeft,CblasLower,CblasTrans,CblasUnit,N,M,a,A1, N, B, N);

                // backward Ux = y
                //cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,N,M,a,A1, N, B, N);

        }else {
                printf("Please input the size of the MATRIX.\n");
        }
}
