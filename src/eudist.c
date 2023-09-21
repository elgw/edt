#include "eudist.h"

static void pass12(const double * restrict B, double * restrict D, const size_t L, const double dx)
/* Pass 1 and 2 for a line (stride 1, since only run along the first dimension) */
{

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
    return;
}


static void pass34(double * restrict D, // Read distances
                   double * restrict D0, // Write distances
                   int * restrict S, int * restrict T, // Temporary pre-allocated storage
                   const int L, // Number of elements in this dimension
                   const int stride,
                   const double d) // voxel size in this dimension
{

    // Make a copy of D into D0
    for(int kk = 0; kk<L; kk++)
    {
        D0[kk] = D[kk*stride];
    }

    // 3: Forward
    int q = 0;
    double w = 0;
    S[0] = 0;
    T[0] = 0;
    const double d2 = d*d;

    for(int u = 1; u<L; u++) // For each column
    {
        // f(t[q],s[q]) > f(t[q], u)
        // f(x,i) = (x-i)^2 + g(i)^2
        while(q >= 0 && ( (pow(d*(T[q]-S[q]), 2) + pow(D0[S[q]], 2)) >=
                          (pow(d*(T[q]-u), 2) +    pow(D0[u], 2)) ) )
            q--;

        if(q<0)
        {
            q = 0;
            S[0] = u;
            T[0] = 0;
        }
        else
        {
            // w = 1 + Sep(s[q],u)
            // Sep(i,u) = (u^2-i^2 +g(u)^2-g(i)^2) div (2(u-i))
            // where division is rounded off towards zero

            w = 1 + trunc(
                          ( pow(d*(double) u,2)  - pow(d*(double) S[q],2)
                            + pow(D0[u],2) - pow(D0[S[q]],2) )
                          /(d2*2.0*(double )(u - S[q]))
                          );

            if(w<L)
            {
                q++;
                S[q] = (int) u; // The point where the segment is a minimizer
                T[q] = (int) w; // The first pixel of the segment
            }
        }
    }

    // 4: Backward
    for(int u = L-1; u > -1 ; u--)
    {
        D[u*stride] = sqrt(pow(d*(u-S[q]), 2) + pow(D0[S[q]], 2));
        if(u <= (int) T[q])
        { q--; }
    }
    return;
}


void
edt(const double * restrict B,
    double * restrict D,
    const size_t M, const size_t N, const size_t P,
    const double dx, const double dy, const double dz)
{

    /* Size of buffers */
    size_t nL = M;
    N > nL ? nL = N : 0;
    P > nL ? nL = P : 0;

#pragma omp parallel
    {
        int * S = calloc(nL, sizeof(int));
        int * T = calloc(nL, sizeof(int));
        double * D0 = calloc(nL, sizeof(double));

#pragma omp for
        for(size_t kk = 0; kk < N*P; kk++) // For each column
        {
            size_t offset = kk*M;
            pass12(B+offset, D+offset, M, dx);
        }


        for(size_t kk = 0; kk < P; kk++) // slice
        {
#pragma omp for
            for(size_t ll = 0; ll<M; ll++) // row
            {
                size_t offset = kk*M*N + ll;
                pass34(D+offset, D0, S, T, N, M, dy);
            }
        }

        if(P > 1)
        {
        for(size_t kk = 0; kk<M; kk++)
        {
#pragma omp for
            for(size_t ll = 0; ll<N; ll++)
            {
                size_t offset = kk + ll*M;
                pass34(D+offset, D0, S, T, P, M*N, dz);
            }
        }
        }

        free(T);
        free(S);
        free(D0);
    }

    return;
}
