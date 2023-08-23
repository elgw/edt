#ifndef eudist_h_
#define eudist_h_

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <stdint.h>
#include "eudist.h"

#define EDT_VERSION_MAJOR "0"
#define EDT_VERSION_MINOR "1"
#define EDT_VERSION_PATCH "0"
#define edt_version EDT_VERSION_MAJOR "."  \
    EDT_VERSION_MINOR "."                        \
    EDT_VERSION_PATCH


/*  Euclidean distance transform
    B specifies a binary mask, 1 == object, 0 = background
    Distances are stored in D
    Matrices are of size M x N x P
    nThreads has to be at least 1.
*/

void edt(double * restrict B,
    double * restrict D,
    const size_t M, const size_t N, const size_t P,
    const double dx, const double dy, const double dz,
    int nThreads);

#endif
