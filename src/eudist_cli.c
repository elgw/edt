#include "eudist.h"

static void edt_brute_force(double * B, double * D, // binary mask and distance
                            size_t M, size_t N, size_t P, // domain size
                            double dx, double dy, double dz) // voxel size
/* Brute force O(n^2) implementation of the euclidean distance
 * transform
 *
 * This is used for testing correctness
 *
 * */
{
    for(size_t mm = 0; mm<M; mm++)
        for(size_t nn = 0; nn<N; nn++)
            for(size_t pp = 0; pp<P; pp++)
            {
                double d_min = INFINITY; // Fallback, everything will be set to this if B is all zeros
                size_t elm  =mm+ M*nn +M*N*pp; // element to consider
                if(B[elm] == 1)
                {
                    D[elm] = 0;
                }
                else
                {
                    for(size_t kk = 0; kk<M; kk++)
                        for(size_t ll = 0; ll<N; ll++)
                            for(size_t qq = 0; qq<P; qq++)
                            {
                                if(B[kk+ M*ll+M*N*qq] == 1)
                                {
                                    double d = pow(dx*((double) kk-(double) mm),2)+
                                        pow(dy*((double) ll-(double) nn),2)+
                                        pow(dz*((double) qq-(double) pp),2);
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


static size_t randi(size_t max)
/* return a value in [0,max] */
{
    size_t val = (size_t) ( (double) rand()/ (double) RAND_MAX * max );
    return val;
}

static double randf(double min, double max)
/* Returns a random value in [min,max] */
{
    return min + (double) rand() / (double) RAND_MAX * (max-min);
}


static void matrix_show(double * M, size_t m, size_t n, size_t p)
/* Utility to print out small matrices */
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
    return;
}

static int test_size(size_t M, size_t N, size_t P, double dx, double dy, double dz, int nThreads)
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
    size_t nB = randi(M*N*P);
    printf("Setting %zu random elements to 1 in B\n", nB);
    for(size_t bb = 0; bb<nB; bb++)
        B[randi(M*N*P-1)] = 1;

    //B[5*2+2] = 1;
    //B[3] = 1;
    if(0)
    {
        printf("Binary mask:\n");
        matrix_show(B, M, N, P);
    }

    //  printf("Edt^2:\n");
    // matrix_show(D, M, N, P);
    // printf("Edt^2 -- brute force reference:\n");
    clock_gettime(CLOCK_MONOTONIC, &start1);
    int bf_run = 0;
    if(M*N*P<1000000)
    {
        bf_run = 1;
        edt_brute_force(B, D_bf,
                        M, N, P,
                        dx, dy, dz);
    }

    clock_gettime(CLOCK_MONOTONIC, &end1);

    clock_gettime(CLOCK_MONOTONIC, &start0);
    edt(B, D,
        M, N, P,
        dx, dy, dz, nThreads);

    clock_gettime(CLOCK_MONOTONIC, &end0);

    double elapsed0, elapsed1;
    elapsed0 = (end0.tv_sec - start0.tv_sec);
    elapsed0 += (end0.tv_nsec - start0.tv_nsec) / 1000000000.0;
    elapsed1 = (end1.tv_sec - start1.tv_sec);
    elapsed1 += (end1.tv_nsec - start1.tv_nsec) / 1000000000.0;

    if(0)
    {
        printf("D_bf:\n");
        matrix_show(D_bf, M, N, P);
    }

    int failed = 0;
    double max_error = 0;

    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        double err = fabs(D[kk] - D_bf[kk]);
        if(err > max_error)
            max_error = err;

        if(err > 10e-5)
            failed = 1;
    }

    if(bf_run)
    {
        if(failed)
        {
            printf("Wrong result! ");
            printf(" -- Largest difference: %f\n", max_error);
        } else
        {
            printf("Correct! ");
            printf("Timing: Edt: %f s, bf: %f s\n", elapsed0, elapsed1);
        }
    }
    else {
        printf("Not verified against brute force\n");
        printf("Timing: %f s\n", elapsed0);
    }


    free(D_bf);
    free(D);
    free(B);

    printf("\n");
    fflush(stdout);

    return failed;
}


static void usage()
{
    printf("Usage:\n");
    printf("-n nThreads : Specify the number of threads to use\n");
    printf("-M #        : Size along first dimension\n");
    printf("-N #        : Size along second dimension\n");
    printf("-P #        : Size along third dimension\n");
    printf("-x #        : Pixel size in first dimension\n");
    printf("-y #        : Pixel size in second dimension\n");
    printf("-z #        : Pixel size in third dimension\n");
    printf("-r          : Run tests with random image and pixel size\n");
    printf("-R          : Use a random seed\n");
    printf("-h          : Show this help message\n");
    return;
}


int main(int argc, char ** argv)
{

    // Defaults:
    int nThreads = 0; // Auto
    size_t M = 1024; size_t N = 1024; size_t P = 60;
    double dx = 1; double dy = 1; double dz = 1;
    int test_one = 1;

    if(argc == 1)
    {
        usage();
        return 0;
    }
    char ch;
    while((ch = getopt(argc, argv, "Rn:M:N:P:rx:y:z:h\n")) != -1)
    {
        switch(ch) {
        case 'n':
            nThreads = atoi(optarg);
            break;
        case 'M':
            M = atoi(optarg);
            break;
        case 'N':
            N = atoi(optarg);
            break;
        case 'P':
            P = atoi(optarg);
            break;
        case 'r':
            test_one = 0;
            break;
        case 'x':
            dx = atof(optarg);
            break;
        case 'y':
            dy = atof(optarg);
            break;
        case 'z':
            dz = atof(optarg);
            break;
        case 'h':
            usage();
            break;
        case 'R':
            srand((unsigned int) time(NULL));
            break;

        }
    }

    if(nThreads < 1)
    {
        printf("Using automatic number of threads\n");
    } else {
        printf("Using %d threads\n", nThreads);
    }

    if(test_one == 1)
    {
        test_size(M,N,P, dx, dy, dz, nThreads);
        return 0;
    }

    printf("Testing random sizes and voxel sizes\n");
    printf("Abort with Ctrl+C\n");

    size_t nTest = 0;
    while(1)
    {
        nTest++;
        printf(" --> Test %zu\n", nTest);
        M = randi(15)+5;
        N = randi(15)+5;
        P = randi(45)+5;
        dx = randf(0.1, 20); // wrong result when dx != 1
        dy = randf(0.1, 20);
        dz = randf(0.1, 20);

        if(test_size(M,N,P, dx, dy, dz, nThreads) > 0)
        {
            printf(" --! Test %zu failed\n", nTest);
            printf("Wrong result for M=%zu, N=%zu, P=%zu, dx=%f, dy=%f, dz=%f\n",
                   M, N, P, dx, dy, dz);
            assert(0);
        }
    }

    return 0;
}
