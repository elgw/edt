#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// #include <omp.h>

#ifndef verbose
#define verbose 2
#endif

size_t randi(size_t max)
{ // return a value in [0,max]

  size_t val = (size_t) ( (double) rand()/ (double) RAND_MAX * max );

  return val;
}

double randf(double min, double max)
{
  return min + (double) rand() / (double) RAND_MAX * (max-min);
}


void edt_brute_force(double * B, double * D, size_t M, size_t N, size_t P, double dx, double dy, double dz)
{

  /* Brute force O(n^2) implementation 
   *
   * This is used for testing correctness
   *
   * */

  for(int mm = 0; mm<M; mm++)
    for(int nn = 0; nn<N; nn++)
      for(int pp = 0; pp<P; pp++)
      {
        double d_min = INFINITY; // Fallback, everything will be set to this if B is all zeros
        size_t elm  =mm+ M*nn +M*N*pp; // element to consider
        if(B[elm] == 1)
        {
          D[elm] = 0;
        }
        else
        {
          for(int kk = 0; kk<M; kk++)
            for(int ll = 0; ll<N; ll++)
              for(int qq = 0; qq<P; qq++)
              {
                if(B[kk+ M*ll+M*N*qq] == 1)
                {
                  double d = pow(dx*(kk-mm),2)+pow(dy*(ll-nn),2)+pow(dz*(qq-pp),2);
                  if( d < d_min)
                    if( ! ((kk == mm ) && (ll == nn) && (qq == pp)))
                    {
                      d_min = d;
                    }
                }
              }
          D[elm] = sqrt(d_min);
        }
      }
return;
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

void matrix_show(double * M, size_t m, size_t n, size_t p)
{
  if(m>20 || n>20 || p>20)
  {
    printf("%zu x %zu x %zu matrix (not shown)\n", m, n, p);
    return;
  }

  for(size_t zz = 0; zz<p; zz++)    
  {
    printf("z=%zu\n", zz);
    for(size_t kk = 0; kk<m ; kk++)
    {
      for(size_t ll = 0; ll<n ; ll++)
      {
        size_t pos = ll*m+kk + zz*m*n;
        printf("%3.2f ", M[pos]);
      }
      printf("\n");
    }
  }
}


void pass12(double * restrict B, double * restrict D, const size_t L, const double dx)
{
  /* Pass 1 and 2 for a line (stride 1, since only run along the first dimension) */

  // Pass 1
  double d = INFINITY;
  for( size_t ll = 0; ll<L; ll++) // For each row
  {
    d = d + dx;
    if(B[ll] == 1)
      d = 0;
    D[ll] = d;
  }

  // Pass 2
  d = INFINITY;
  for( size_t ll = L ; ll-- > 0; ) // For each row
  {
    d = d + dx;
    if(B[ll] == 1)
      d = 0;
    if(d<D[ll])
      D[ll] = d;
  }

}

void pass34(double * restrict D, double * restrict D0, 
    int * restrict S, double * restrict T, 
    const int L, const int stride, const double d)
{
  // 3: Forward
  int q = 0;
  double w = 0;
  S[0] = 0;
  T[0] = 0;
  double d2 = d*d;

  for(int u = 1; u<L; u++) // For each column
  {
    // f(t[q],s[q]) > f(t[q], u)
    // f(x,i) = (x-i)^2 + g(i)^2
    while(q >= 0 && ( (pow(d*(T[q]-S[q]),2) + pow(D[stride*S[q]],2)) > 
                      (pow(d*(T[q]-u), 2) +  pow(D[stride*u],2)) ) )
      q--;

    if(q<0)
    {
      if(verbose > 2)
        printf("reset\n");
      q = 0;
      S[0] = u;
    }
    else
    {
      // w = 1 + Sep(s[q],u)
      // Sep(i,u) = (u^2-i^2 +g(u)^2-g(i)^2) div (2(u-i))
      // where division is rounded off towards zero
      //
      // TODO: derive correctly for non-unit step length
      w = 1 + floor0( ( pow(d*u,2)  - pow(d*S[q],2) 
                     + pow(D[stride*u],2) - pow(D[stride*S[q]],2))/(2*d2*(u-S[q])));
      // because of overflow, w is double. T does not have to be
      // double
      if(verbose > 3)
        printf("u/kk: %d, S[q] = %d, q: %d w: %f\n", u, S[q], q, w);

      if(w<L)
      {
        q++;
        S[q] = (int) u;
        T[q] = (int) w;
      }
    }

  }

  // 4: Backward  
  for(int u = L-1; u > -1 ; u--)
  {
    //dt[u,y]:=f(u,s[q])
    if(verbose>3)
      printf("u: %d, q: %d S[%d] = %d\n", u, q, q, S[q]);
    D[u*stride] = sqrt(pow(d*(u-S[q]),2)+pow(D0[stride*S[q]], 2));
    if(u == T[q])
      q--;
  }
}


void edt(double * B, double * D, size_t M, size_t N, size_t P, 
    double dx, double dy, double dz)
{
  // Euclidean distance transform from objects in M, into D
  // Matrices are of size M x N

  size_t nL = max_size_t(M,max_size_t(N, P));

  // 
  // First dimension, pass 1 and 2
  //
//#pragma omp parallel for
  for(size_t kk = 0; kk<N*P; kk++) // For each column
  {
    size_t offset = kk*M;
    pass12(B+offset, D+offset, M, dx);
  }

  if(verbose>1)
  {
    printf("D:\n");
    matrix_show(D, M, N, P);
  }

  int * S = malloc(nL*nL*sizeof(double));
  double * T = malloc(nL*nL*sizeof(double));

  //maintain only one line, not all of this
  double * D0 = malloc(M*N*P*sizeof(double));
  memcpy(D0, D, M*N*P*sizeof(double));

  // 
  // Second dimension, Pass 3 and 4
  //

  int length = N;
  int stride = M;

//#pragma omp parallel for
  for(int kk = 0; kk<P; kk++) // slice
  {
    for(int ll = 0; ll<M; ll++) // row
    {
      size_t offset = kk*M*N + ll;
      // TODO: does the threads get separate copies of S and T?
      pass34(D+offset, D0+offset, S+nL*kk, T+nL*kk, length, stride, dy);
    }
  }

  if(verbose>1)
  {
    printf("D2:\n");
    matrix_show(D, M, N, P);
  }


  // Third dimension
  memcpy(D0, D, M*N*P*sizeof(double)); // TODO: Threaded

  if(P>1)
  {
    length = P;
    stride = M*N;
//#pragma omp parallel for
    for(int kk = 0; kk<M; kk++)
    {
      for(int ll = 0; ll<N; ll++)
      {
        //        int kk = 2; int ll = 0;
        size_t offset = kk + ll*M;
        //printf("O: %zu S: %zu L: %zu\n", offset, stride, length);
        pass34(D+offset, D0+offset, S+nL*kk, T+nL*kk, length, stride, dz);
      }
    }
  }

  if(verbose>1)
  {
    printf("D3:\n");
    matrix_show(D, M, N, P);
  }


  free(D0);
  free(T);
  free(S);

}


int test_size(size_t M, size_t N, size_t P, double dx, double dy, double dz)
{
  printf("Problem size: %zu x %zu x %zu\n", M, N, P);
  printf("Voxel size: %.2f x %.2f x %.2f\n", dx, dy, dz);

  /* Allocate memory */
  double * B = calloc(M*N*P, sizeof(double));
  double * D = calloc(M*N*P, sizeof(double));
  double * D_bf = calloc(M*N*P, sizeof(double));

  // For timing
  struct timespec start0, end0, start1, end1;

  /* Initialize binary mask */
  size_t nB = randi(10);
  printf("Setting %zu random elements to 1 in B\n", nB);
  for(int bb = 0; bb<nB; bb++)
    B[randi(M*N*P)-1] = 1;

  //B[5*2+2] = 1;
  //B[3] = 1;
  printf("Binary mask:\n");
  matrix_show(B, M, N, P);

  //  printf("Edt^2:\n");
  // matrix_show(D, M, N, P);
  // printf("Edt^2 -- brute force reference:\n");
  clock_gettime(CLOCK_MONOTONIC, &start1);
  edt_brute_force(B, D_bf, 
      M, N, P, 
      dx, dy, dz);

  clock_gettime(CLOCK_MONOTONIC, &end1);

  clock_gettime(CLOCK_MONOTONIC, &start0);
  edt(B, D, 
      M, N, P,
      dx, dy, dz);

  clock_gettime(CLOCK_MONOTONIC, &end0);

  double elapsed0, elapsed1;
  elapsed0 = (end0.tv_sec - start0.tv_sec);
  elapsed0 += (end0.tv_nsec - start0.tv_nsec) / 1000000000.0;
  elapsed1 = (end1.tv_sec - start1.tv_sec);
  elapsed1 += (end1.tv_nsec - start1.tv_nsec) / 1000000000.0;

  printf("D_bf:\n");
  matrix_show(D_bf, M, N, P);

  int wrong_result = 0;
  double max_error = 0;

  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    double err = fabs(D[kk] - D_bf[kk]); 
    if(err > max_error)
      max_error = err;
      
    if(max_error > 10e-5)
    {
      wrong_result = 1;
    }
  }

    printf("Largest difference: %f\n", max_error);

  if(wrong_result)
  {
    printf("Wrong result! ");
  } else
  {
    printf("Correct! ");
  }

  printf("Edt: %f s, bf: %f s\n", elapsed0, elapsed1);

  free(D_bf);
  free(D);
  free(B);

  return wrong_result;
}

int main(int argc, char ** argv)
{

  printf("Testing random sizes and voxel sizes\n");
  printf("TODO: Present a summary when aborting\n");
  printf("Abort with Ctrl+C\n");

//  test_size(8, 9, 1, 1.812511, 3.9250, 13.298);
//  test_size(8, 9, 1, 1, 2, 1);

//  return 0;

  size_t nTest = 0;
  while(1)
  {
    nTest++;
    printf(" --> Test %zu\n", nTest);
    size_t M = 3; //randi(15)+5;
    size_t N = randi(15)+5;
    size_t P = 1; // randi(45)+5;
    double dx = 2.1;//randf(0.1, 20); // wrong result when dx != 1
    double dy = 1;//randf(0.1, 20);
    double dz = 1;//randf(0.1, 20);
    
    if(test_size(M,N,P, dx, dy, dz) > 0)
    {
    printf(" --! Test %zu failed\n", nTest);
      printf("Wrong result for M=%zu, N=%zu, P=%zu, dx=%f, dy=%f, dz=%f\n",
          M, N, P, dx, dy, dz);
      assert(0);
    }
  }

  return 0;
}
