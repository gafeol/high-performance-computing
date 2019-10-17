/* ==================================================================== *
 *								        *
 *  block_dgemm.c -- Implemant a block matrix multiplication routine    *
 *                                                                      *
 * ==================================================================== */

#include "square_dgemm.h"

#include <stdio.h>
#include <stdlib.h>

/* block parameter ... */
#ifndef BLOCK_SIZE
#  define BLOCK_SIZE ((unsigned) 16)
#endif


/**
 *  square_dgemm -- multiply two block matrices A and B adding result to C, result is C = C + A*B
 */

void slow_read(double* to, const double *from, int i, int j, const unsigned M){
    int cnt = 0;
    for(int ii = i;ii < i + BLOCK_SIZE;++ii){
        for(int jj = j;jj < j + BLOCK_SIZE;++jj){
            if(ii >= M || jj >= M)
                to[cnt++] = 0;
            else
                to[cnt++] = from[ii * M + jj];
        }
    }
}

void naive_mm(double *C, double *A, double *B){
    for(int i=0;i<BLOCK_SIZE;++i){
        for(int j=0;j<BLOCK_SIZE;++j){
            double ans = 0;
            for(int k=0;k<BLOCK_SIZE;++k)
                ans += A[i*BLOCK_SIZE + k] * B[k*BLOCK_SIZE + j];
            C[i * BLOCK_SIZE + j] += ans;
        }
    }
}

void slow_write(double *C, double *Cb, int i, int j, const unsigned M){
    for(int ii=0;ii < BLOCK_SIZE && i + ii < M;++ii){
        for(int jj=0;jj < BLOCK_SIZE && j + jj < M;++jj){
            C[(i+ii) * M + j + jj] += Cb[ii*BLOCK_SIZE + jj];
        }
    }
}

void square_dgemm (const double  *A, const double  *B,  double  *C, const unsigned  M){
    double Ab[BLOCK_SIZE * BLOCK_SIZE], Bb[BLOCK_SIZE * BLOCK_SIZE], Cb[BLOCK_SIZE * BLOCK_SIZE];
    for (int i=0;i<M; i+=BLOCK_SIZE) {
        for (int j=0;j<M; j+=BLOCK_SIZE) {
            slow_read(Cb, C, i, j, M);
            for(int k=0;k<M;k+=BLOCK_SIZE){
                slow_read(Ab, A, i, k, M);
                slow_read(Bb, B, k, j, M);
                naive_mm(Cb, Ab, Bb);
            }
            slow_write(C, Cb, i, j, M);
        }
    }
}
