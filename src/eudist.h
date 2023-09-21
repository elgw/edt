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
#define EDT_VERSION_PATCH "1"
#define edt_version EDT_VERSION_MAJOR "."       \
    EDT_VERSION_MINOR "."                       \
    EDT_VERSION_PATCH


/* Euclidean distance transform implemented from
 * https://doi.org/10.1007/0-306-47025-X_36
 *
 * Input:
 *    B specifies a binary mask, 1 == object, 0 = background
 *    Matrices are of size M x N x P
 *    dx, dy, dz specifies the pixel size
 *    If nThreads < 1 the number of threads will be set
 *    as sysconf(_SC_NPROCESSORS_ONLN)/2
 *
 * Output:
 *    Distances are stored in D
*/

void edt(const double * restrict B,
         double * restrict D,
         const size_t M, const size_t N, const size_t P,
         const double dx, const double dy, const double dz,
         int nThreads);

#endif
