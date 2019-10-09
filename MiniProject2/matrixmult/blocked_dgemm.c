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
#  define BLOCK_SIZE ((unsigned) 1)
#endif



/**
 *  square_dgemm -- multiply two block matrices A and B adding result to C, result is C = C + A*B
 */
void square_dgemm (const double  *A, const double  *B,  double  *C, const unsigned  M)
{
    for (int i = 0; i < M; i+=BLOCK_SIZE) {
        for (int j = 0; j < M; j+=BLOCK_SIZE) {
            for(int k=0;k<M;k+=BLOCK_SIZE){
                for(int ii=i;ii<i+BLOCK_SIZE;ii++){
                    for(int jj=j;jj<j+BLOCK_SIZE;jj++){
                        C[ii*M + jj] += A[ii*M + k] * B[k*M + jj];
                    }
                }
            }
        }
    }
}
