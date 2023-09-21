#include "eudist.h"

/* The structure for the "private" given to each thread */
typedef struct{
    const double * B; // Binary mask
    double * D; // Output
    double * D0; // Line buffer of size max(M, N, P)

    int * S; // Line buffers for pass 3 and 4
    int * T;

    int thrId;
    int nThreads;

    size_t M; // Image dimensions
    size_t N;
    size_t P;

    double dx; // Voxel size
    double dy;
    double dz;
} thrJob;

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


static void * pass12_t(void * data)
{
    thrJob * job = (thrJob *) data;

    int last = job->N*job->P;
    int from = job->thrId * last/job->nThreads;
    int to = (job->thrId + 1)*last/job->nThreads-1;

//  printf("Thread: %d. From: %d To: %d, last: %d\n", job->thrId, from, to, last);
//  fflush(stdout);

    for(int kk = from; kk<=to; kk++) // For each column
    {
        size_t offset = kk*job->M;
        pass12(job->B+offset, job->D+offset, job->M, job->dx);
    }
    return NULL;
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

            w = 1 + trunc( ( pow(d*u,2)  - pow(d*(double) S[q],2)
                             + pow(D0[u],2) - pow(D0[S[q]],2))/(d2*2.0*(u-(double) S[q]))
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
        //dt[u,y]:=f(u,s[q])

        D[u*stride] = sqrt(pow(d*(u-S[q]), 2) + pow(D0[S[q]], 2));
        if(u <= (int) T[q])
        { q--; }
    }
    return;
}

static void * pass34y_t(void * data)
{
    thrJob * job = (thrJob *) data;

    // Second dimension
    int length = job->N;
    int stride = job->M;
    double dy = job->dy;

    for(size_t kk = 0; kk< job->P; kk++) // slice
    {
        for(size_t ll = job->thrId; ll<job->M; ll=ll+job->nThreads) // row
        {
            size_t offset = kk*job->M*job->N + ll;
            pass34(job->D+offset, job->D0, job->S, job->T, length, stride, dy);
        }
    }
    return NULL;
}

static void * pass34z_t(void * data)
{
    thrJob * job = (thrJob *) data;

    // Second dimension
    int length = job->P;
    int stride = job->M*job->N;
    double dz = job->dz;

    for(size_t kk = job->thrId; kk<job->M; kk=kk+job->nThreads)
    {
        for(size_t ll = 0; ll<job->N; ll++)
        {
            size_t offset = kk + ll*job->M;
            pass34(job->D+offset, job->D0, job->S, job->T, length, stride, dz);
        }
    }
    return NULL;
}


void
edt(const double * restrict B,
    double * restrict D,
    const size_t M, const size_t N, const size_t P,
    const double dx, const double dy, const double dz,
    int nThreads)
{

    size_t nL = M;
    N > nL ? nL = N : 0;
    P > nL ? nL = P : 0;

    if(nThreads < 1)
    {
        nThreads = sysconf(_SC_NPROCESSORS_ONLN)/2;
    }

#ifdef timings
    struct timespec tic, toc;
    double tot;
#endif

    /* Set up threads and their buffers */

    pthread_t threads[nThreads];
    thrJob jobs[nThreads];

    for(int kk = 0; kk<nThreads; kk++)
    {
        jobs[kk].B = B;
        jobs[kk].D = D;
        jobs[kk].S = calloc(nL, sizeof(int));
        jobs[kk].T = calloc(nL, sizeof(int));
        jobs[kk].D0 = calloc(nL, sizeof(double));
        jobs[kk].thrId = kk;
        jobs[kk].nThreads = nThreads;
        jobs[kk].M = M;
        jobs[kk].N = N;
        jobs[kk].P = P;
        jobs[kk].dx = dx;
        jobs[kk].dy = dy;
        jobs[kk].dz = dz;
    }


    /* First dimension, pass 1 and 2 */

#ifdef timings
    clock_gettime(CLOCK_MONOTONIC, &tic);
#endif

    for(int kk = 0; kk<nThreads; kk++) {
        pthread_create(&threads[kk], NULL, pass12_t, &jobs[kk]); }

    for(int kk = 0; kk<nThreads; kk++) {
        pthread_join(threads[kk], NULL); }


#ifdef timings
    clock_gettime(CLOCK_MONOTONIC, &toc);
    tot = (toc.tv_sec - tic.tv_sec);
    tot += (toc.tv_nsec - tic.tv_nsec) / 1000000000.0;
    printf("x Took %f s\n", tot);
#endif


    /* Second dimension */
#ifdef timings
    clock_gettime(CLOCK_MONOTONIC, &tic);
#endif

    for(int kk = 0; kk<nThreads; kk++) {
        pthread_create(&threads[kk], NULL, pass34y_t, &jobs[kk]); }

    for(int kk = 0; kk<nThreads; kk++) {
        pthread_join(threads[kk], NULL); }

#ifdef timings
    clock_gettime(CLOCK_MONOTONIC, &toc);
    tot = (toc.tv_sec - tic.tv_sec);
    tot += (toc.tv_nsec - tic.tv_nsec) / 1000000000.0;
    printf("y Took %f s\n", tot);
#endif

    /* Third dimension */
#ifdef timings
    clock_gettime(CLOCK_MONOTONIC, &tic);
#endif
    if(P>1)
    {
        for(int kk = 0; kk<nThreads; kk++) {
            pthread_create(&threads[kk], NULL, pass34z_t, &jobs[kk]); }

        for(int kk = 0; kk<nThreads; kk++) {
            pthread_join(threads[kk], NULL); }
    }
#ifdef timings
    clock_gettime(CLOCK_MONOTONIC, &toc);
    tot = (toc.tv_sec - tic.tv_sec);
    tot += (toc.tv_nsec - tic.tv_nsec) / 1000000000.0;
    printf("z Took %f s\n", tot);
#endif

    /* Finalize */
    for(int kk = 0; kk<nThreads; kk++)
    {
        free(jobs[kk].S);
        free(jobs[kk].T);
        free(jobs[kk].D0);
    }

    return;
}
