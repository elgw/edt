#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


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


void edt2(double * B, double * D, size_t M, size_t N)
{
  // Euclidean distance transform from objects in M, into D
  // Matrices are of size M x N

  size_t nL = max_size_t(M,N);
  int bignum = 2*nL+1;

  // Pass 1
  for(size_t kk = 0; kk<N; kk++)
  {
    double d = bignum;
    for( size_t ll = 0; ll<M; ll++)
    {
      size_t pos = kk*M+ll;
      d++;
      if(B[pos] == 1)
        d = 0;
      D[pos] = d;
    }
  }

  printf("D1:\n");
  matrix_show(D, M, N);

  // Pass 2
  for(size_t kk = N ; kk-- > 0 ; )
  {
    double d = bignum;
    for( size_t ll = M ; ll-- > 0; )
    {
      size_t pos = kk*M+ll;
      d++;
      if(B[pos] == 1)
        d = 0;
      if(d<D[pos])
        D[pos] = d;
    }
  }


  int * S = malloc(nL*sizeof(double));
  int * T = malloc(nL*sizeof(double));
  printf("D2:\n");
  matrix_show(D, M, N);

  //maintain only one line, not all of this
  double * D0 = malloc(M*N*sizeof(double));
  memcpy(D0, D, M*N*sizeof(double));


  // Pass 3 and 4
  for(int ll = 0; ll<N; ll++)
  {
    // 3: Forward
    int q = 0;
    int w = 0;
    S[0] = 0;
    T[0] = 0;

    for(int kk = 1; kk<M; kk++)
    {
      // f(t[q],s[q]) > f(t[q], u)
      // f(x,i) = (x-i)^2 + g(i)^2
      while(q>=0 && ((pow(T[q]-S[q],2) + pow(D[M*S[q]+ll],2)) > (pow(T[q]-kk, 2) + pow(D[M*kk+ll],2))))
        q--;

      if(q<0)
      {
        printf("reset\n");
        q = 0;
        S[0] = kk;
      }
      else
      {
        // w = 1 + Sep(s[q],u)
        // Sep(i,u) = (u^2-i^2+g(u)^2-g(i)^2) div (2(u-i))
        w = 1 + floor0((pow(kk,2)-pow(S[q],2) 
              + pow(D[M*kk+ll],2) - pow(D[M*S[q]+ll],2))/(2*(kk-S[q])));
        printf("u/kk: %d, S[q] = %d, q: %d w: %d\n", kk, S[q], q, w);

        if(w<(int) M)
        {
          q++;
          S[q] = kk;
          T[q] = w;
        }
      }

    }

    // 4: Backward  
    for(int kk = M-1; kk > -1 ; kk--)
    {
      //dt[u,y]:=f(u,s[q])
      printf("kk: %d, q: %d S[%d] = %d\n", kk, q, q, S[q]);
      D[kk*M+ll] = pow(kk-S[q],2)+pow(D0[M*S[q]+ll], 2);
      if(kk == T[q])
        q--;
    }
  }

  free(D0);
  free(T);
  free(S);

}


int main(int argc, char ** argv)
{

  size_t M = 5;
  size_t N = 5;

  double * B = calloc(M*N, sizeof(double));
  double * D = calloc(M*N, sizeof(double));

  B[5*2+2] = 1;
  B[3] = 1;
  printf("B:\n");
  matrix_show(B, M, N);

  printf("D:\n");
  edt2(B, D, M, N);

  matrix_show(D, M, N);
  return 0;
}
