#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


size_t max_size_t(size_t a, size_t b)
{
  if(a>b)
    return a;
  return b;
}


void edt2(double * B, double * D, size_t M, size_t N)
{
  // Euclidean distance transform from objects in M, into D
  // Matrices are of size M x N


  // Pass 1
  for(size_t kk = 0; kk<N; kk++)
  {
    double d = INFINITY;
    for( size_t ll = 0; ll<M; ll++)
    {
      size_t pos = kk*M+ll;
      d++;
      if(B[pos] == 1)
        d = 0;
      D[pos] = d;
    }
  }

  // Pass 2
  for(size_t kk = N-1 ; kk-->0 ; )
  {
    double d = INFINITY;
    for( size_t ll = M-1 ; ll-- > 0; )
    {
      size_t pos = kk*M+ll;
      d++;
      if(B[pos] == 1)
        d = 0;
      if(d<D[pos])
        D[pos] = d;
    }
  }


  size_t nL = max_size_t(M,N);
  double * S = malloc(nL*sizeof(double));
  double * T = malloc(nL*sizeof(double));


  // Pass 3 and 4
  for(size_t ll = 0; ll<N; ll++)
  {
    // 3: Forward
    int q = 0;
    double w = 0;
    S[0] = 0;
    T[0] = 0;

    for(size_t kk = 0; kk<M; kk++)
    {
      if(q<0)
      {
        q = 0;
        S[0] = kk;
      }
      else
      {
        w = 1 + (pow(ll,2)-pow(S[kk],2)); //... sep(s[kk], u);
      }
      if(w<M)
      {
        q++;
        S[kk] = kk;
        T[q] = w;
      }

      }

    // 4: Backward  
    for(size_t kk = M-1; kk-- > 0 ;)
    {
      D[kk*M+ll] = pow(kk-kk,2)+pow(S[q], 2);
      if(kk == T[q])
        q--;
    }
  }

  free(T);
  free(S);

}

void matrix_show(double * M, size_t m, size_t n)
{
  for(size_t kk = 0; kk<m ; kk++)
  {
    for(size_t ll = 0; ll<n ; ll++)
    {
      size_t pos = ll*m+kk;
      printf("%f ", M[pos]);
    }
    printf("\n");
  }
}


int main(int argc, char ** argv)
{

  size_t M = 5;
  size_t N = 5;

  double * B = calloc(M*N, sizeof(double));
  double * D = calloc(M*N, sizeof(double));

  B[5*2+2] = 1;
  printf("B:\n");
  matrix_show(B, M, N);
  
  printf("D:\n");
  edt2(B, D, M, N);

  matrix_show(D, M, N);
  return 0;
}
