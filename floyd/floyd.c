#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

double ts, tp;

int square(int n, int* restrict g, int* restrict gnew)
{
    int done = 1;
    #pragma omp parallel for shared(g, gnew) reduction(&& : done)
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            int gij = gnew[j*n+i];
            for (int k = 0; k < n; ++k) {
                int gik = g[k*n+i];
                int gkj = g[j*n+k];
                if (gik + gkj < gij) {
                    gij = gik+gkj;
                    done = 0;
                }
            }
            gnew[j*n+i] = gij;
        }
    }
    return done;
}

static inline void infinitize(int n, int* g)
{
    for (int i = 0; i < n*n; ++i)
        if (g[i] == 0)
            g[i] = 10000;
}

static inline void deinfinitize(int n, int* g)
{
    for (int i = 0; i < n*n; ++i)
        if (g[i] == 10000)
            g[i] = 0;
}

int* shortest_paths_p(int n, int* restrict g)
{
    int* restrict gg = (int*) calloc(n*n, sizeof(int));
    memcpy(gg, g, n*n * sizeof(int));
    tp = omp_get_wtime();
    infinitize(n, gg);
    for (int i = 0; i < n*n; i += n+1)
        gg[i] = 0;
    int* restrict gnew = (int*) calloc(n*n, sizeof(int));
    memcpy(gnew, gg, n*n * sizeof(int));
    for (int done = 0; !done; ) {
        done = square(n, gg, gnew);
        memcpy(g, gnew, n*n * sizeof(int));
    }
    free(gnew);
    deinfinitize(n, g);
    tp = omp_get_wtime() - tp;
    return gg;
}

int* shortest_paths_s(int n, int* restrict g)
{
    int* restrict gnew = (int*) calloc(n*n, sizeof(int));
    memcpy(gnew, g, n*n * sizeof(int));
    ts = omp_get_wtime();
    infinitize(n, gnew);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++) {
                if (gnew[j*n+k] > gnew[j*n+i] + gnew[i*n+k])
                    gnew[j*n+k] = gnew[j*n+i] + gnew[i*n+k];
            }
    }
    for (int i = 0; i < n*n; i += n+1)
        gnew[i] = 0;
    ts = omp_get_wtime() - ts;
    deinfinitize(n, gnew);
    return gnew;
}

int* gen_graph(int n)
{
    int* g = calloc(n*n, sizeof(int));
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            if(j > i)
                g[j*n+i] = (rand() % 100);
            else g[j*n+i] = -1;
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            if(i > j)
                g[j*n+i] = g[i*n+j];
        g[j*n+j] = 0;
    }
    return g;
}

void write_matrix(const char* fname, int n, int* a)
{
    FILE* fp = fopen(fname, "w+");
    if (fp == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", fname);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            fprintf(fp, "%d ", a[j*n+i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int main(int argc, char** argv)
{
    int n = (argc > 1) ? atoi(argv[1]) : 100;
    const char *genname = (argc > 2) ? argv[2] : NULL;
    const char *sername = (argc > 3) ? argv[3] : NULL;
    const char *parname = (argc > 4) ? argv[4] : NULL;
    if (n < 1) {
        fprintf(stderr, "Invalid size of graph");
        exit(EXIT_FAILURE);
    }

    int* g = gen_graph(n);
    if (genname)
        write_matrix(genname,  n, g);

    int* gs = shortest_paths_s(n, g);
    int* gp = shortest_paths_p(n, g);

    printf("OpenMP with %d threads\n", omp_get_max_threads());
    printf("n:     %d\n", n);
    printf("Time ser (sec): %.6f\nTime parallel (sec): %.6f\nS(n): %.6f\n",
                                          ts, tp, ts / tp);

    if (sername)
        write_matrix(sername, n, gs);
    if (parname)
        write_matrix(parname, n, gs);

    free(g);
    free(gs);
    free(gp);
    return 0;
}
