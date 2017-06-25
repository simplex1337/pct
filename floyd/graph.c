#include "graph.h"

graph *graph_create(int nvertices)
{
    struct graph *g;
    g = (struct graph *) malloc(sizeof(*g));
    g->nvertices = nvertices;
    g->m = (int **) malloc(nvertices * sizeof(int *));
    if (g->m != NULL) {
        for (int i = 0; i < nvertices; i++) {
            g->m[i] = (int *) malloc(nvertices * sizeof(int));
        }
    } else {
        printf("Can't create buffer!\n");
    }
    graph_clear(g);             // Опционально, O(n^2)
    return g;
}

void graph_clear(graph * g)
{
    int i, j;
    for (i = 0; i < g->nvertices; i++)
        for (j = 0; j < g->nvertices; j++)
            g->m[i][j] = 100000;
    for (i = 0; i < g->nvertices; i++)
        g->m[i][i] = 0;
}

void graph_free(graph * g)
{
    for (int i = 0; i < g->nvertices; i++)
        free(g->m[i]);
    free(g->m);
    free(g);
}

void graph_set_edge(graph * g, int i, int j, int w)
{
    g->m[i - 1][j - 1] = w;
}

int graph_get_edge(graph * g, int i, int j)
{
    return g->m[i - 1][j - 1];
}

void printm(int **matrix, int cnt)      //вывод матрицы на экран
{
    for (int i = 0; i < cnt; i++) {
        for (int j = 0; j < cnt; j++) {
            if (matrix[i][j] == 100000)
                printf("inf ");
            else
                printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double graph_short_path_floyd_serial(graph * g)
{
    // printm(g->m, g->nvertices);
    double t = omp_get_wtime();
    for (int i = 0; i < g->nvertices; i++) {
        for (int j = 0; j < g->nvertices; j++)
            for (int k = 0; k < g->nvertices; k++) {
                if (g->m[j][k] > g->m[j][i] + g->m[i][k])
                    g->m[j][k] = g->m[j][i] + g->m[i][k];
            }
        //printm(g->m, g->nvertices);
    }
    t = omp_get_wtime() - t;
    return t;
    //return g;
}

double graph_short_path_floyd_parallel(graph * g)
{
    double t = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < g->nvertices; i++) {
        for (int j = 0; j < g->nvertices; j++)
            for (int k = 0; k < g->nvertices; k++) {
                if (g->m[j][k] > g->m[j][i] + g->m[i][k])
                    g->m[j][k] = g->m[j][i] + g->m[i][k];
            }
    }
    t = omp_get_wtime() - t;
    return t;
    //return g;
}
