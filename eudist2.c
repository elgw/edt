#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#define verbose 1

void edt2_brute_force(double * B, double * D, size_t M, size_t N)
{

  /* Brute force O(n^2) implementation 
   *
   * This is used for testing correctness
   *
   * */

  for(int mm = 0; mm<M; mm++)
    for(int nn = 0; nn<N; nn++)
    {
      double d_min = INFINITY; // Fallback, everything will be set to this if B is all zeros
      if(B[mm+M*nn] == 1)
      {
        D[mm+M*nn] = 0;
      }
      else
      {
        for(int kk = 0; kk<M; kk++)
          for(int ll = 0; ll<N; ll++)
          {
            if(B[kk+M*ll] == 1)
            {
              double d = pow(kk-mm,2)+pow(ll-nn,2);
              if( d < d_min)
                if( ! ((kk == mm ) && (ll == nn)) )
                {
                  d_min = d;
                }
            }
          }
        D[mm+M*nn] = d_min;
      }
    }

}

// round towards zero
double floor0( double value )
{
  if (value < 0.0)
    return ceil( value );
  else
    return floor( value );
}

size_t max_size_t(size_t a, size_t b)
{
  if(a>b)
    return a;
  return b;
}

void matrix_show(double * M, size_t m, size_t n)
{
  if(m>20)
    return;
  if(n>20)
    return;

  for(size_t kk = 0; kk<m ; kk++)
  {
    for(size_t ll = 0; ll<n ; ll++)
    {
      size_t pos = ll*m+kk;
      printf("%3.0f ", M[pos]);
    }
    printf("\n");
  }
}


void pass12(double * B, double * D, size_t L)
{
  /* Pass 1 and 2 for a line (stride 1, since only run along the first dimension) */
  
  // Pass 1
  double d = INFINITY;
  for( size_t ll = 0; ll<L; ll++) // For each row
  {
    d++;
    if(B[ll] == 1)
      d = 0;
    D[ll] = d;
  }

  // Pass 2
  d = INFINITY;
  for( size_t ll = L ; ll-- > 0; ) // For each row
  {
    d++;
    if(B[ll] == 1)
      d = 0;
    if(d<D[ll])
      D[ll] = d;
  }

}

void pass34(double * D, double * D0, int * S, double * T, int L, int stride)
{
    // 3: Forward
    int q = 0;
    double w = 0;
    S[0] = 0;
    T[0] = 0;

    for(int kk = 1; kk<L; kk++) // For each column
    {
      // f(t[q],s[q]) > f(t[q], u)
      // f(x,i) = (x-i)^2 + g(i)^2
      while(q>=0 && ((pow(T[q]-S[q],2) + pow(D[stride*S[q]],2)) > (pow(T[q]-kk, 2) + pow(D[stride*kk],2))))
        q--;

      if(q<0)
      {
        if(verbose > 2)
          printf("reset\n");
        q = 0;
        S[0] = kk;
      }
      else
      {
        // w = 1 + Sep(s[q],u)
        // Sep(i,u) = (u^2-i^2+g(u)^2-g(i)^2) div (2(u-i))
        w = 1 + floor0((pow(kk,2)-pow(S[q],2) 
              + pow(D[stride*kk],2) - pow(D[stride*S[q]],2))/(2*(kk-S[q])));
        if(verbose > 2)
          printf("u/kk: %d, S[q] = %d, q: %d w: %f\n", kk, S[q], q, w);

        if(w<(int) L)
        {
          q++;
          S[q] = kk;
          T[q] = w;
        }
      }

    }

    // 4: Backward  
    for(int kk = L-1; kk > -1 ; kk--)
    {
      //dt[u,y]:=f(u,s[q])
      if(verbose>1)
        printf("kk: %d, q: %d S[%d] = %d\n", kk, q, q, S[q]);
      D[kk*stride] = pow(kk-S[q],2)+pow(D0[stride*S[q]], 2);
      if(kk == T[q])
        q--;
    }
}


void edt2(double * B, double * D, size_t M, size_t N)
{
  // Euclidean distance transform from objects in M, into D
  // Matrices are of size M x N

  size_t nL = max_size_t(M,N);

  // 
  // First dimension, pass 1 and 2
  //
  
  for(size_t kk = 0; kk<N; kk++) // For each column
  {
    pass12(B+kk*M, D+kk*M, M);
  }

  if(verbose>1)
  {
    printf("D:\n");
    matrix_show(D, M, N);
  }

  int * S = malloc(nL*sizeof(double));
  double * T = malloc(nL*sizeof(double));

  if(verbose>1)
  {
    printf("D2:\n");
    matrix_show(D, M, N);
  }

  //maintain only one line, not all of this
  double * D0 = malloc(M*N*sizeof(double));
  memcpy(D0, D, M*N*sizeof(double));

  // 
  // Second dimension, Pass 3 and 4
  //
  
  for(int ll = 0; ll<M; ll++) // Select row
  {
    int length = N;
    int stride = M;
    pass34(D+ll, D0+ll, S, T, length, stride);
  }

  free(D0);
  free(T);
  free(S);

}


int main(int argc, char ** argv)
{

  size_t M = 5;
  size_t N = 7;

  printf("Problem size: %zu x %zu\n", M, N);

  double * B = calloc(M*N, sizeof(double));
  double * D = calloc(M*N, sizeof(double));
  double * D_bf = calloc(M*N, sizeof(double));

  // For timing
  struct timespec start0, end0, start1, end1;

  B[5*2+2] = 1;
  B[3] = 1;
  printf("Binary mask:\n");
  matrix_show(B, M, N);

  printf("Edt^2:\n");
  clock_gettime(CLOCK_MONOTONIC, &start0);
  for(int kk = 0 ; kk<1000; kk++)
    edt2(B, D, M, N);
  clock_gettime(CLOCK_MONOTONIC, &end0);

  matrix_show(D, M, N);
  int has_ref = 0;
  printf("Edt^2 -- brute force reference:\n");
  clock_gettime(CLOCK_MONOTONIC, &start1);
  if(M < 100 && N<100)
  {
    edt2_brute_force(B, D_bf, M, N);
    has_ref = 1;
  } else {
    printf("Too large problem, not testing\n");
  }
  clock_gettime(CLOCK_MONOTONIC, &end1);

  double elapsed0, elapsed1;
  elapsed0 = (end0.tv_sec - start0.tv_sec);
  elapsed0 += (end0.tv_nsec - start0.tv_nsec) / 1000000000.0;
  elapsed1 = (end1.tv_sec - start1.tv_sec);
  elapsed1 += (end1.tv_nsec - start1.tv_nsec) / 1000000000.0;

  matrix_show(D_bf, M, N);

  if(has_ref == 1)
  {
    for(size_t kk = 0; kk<M*N; kk++)
      assert(D[kk] == D_bf[kk]);

    printf("Correct result\n");
  }
  else 
  {
    printf("Correctness not tested\n");
  }

  printf("Edt: %f s, bf: %f s\n", elapsed0, elapsed1);

  free(D_bf);
  free(D);
  free(B);

  return 0;
}
