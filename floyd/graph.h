#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

typedef struct graph {
    int nvertices;              /* Число вершин */
    int **m;                    /* Матрица n x n */
} graph;

graph *graph_create(int nvertices);
void graph_clear(graph *g);
void graph_free(graph *g);
void graph_set_edge(graph *g, int i, int j, int w);
int graph_get_edge(graph *g, int i, int j);
double graph_short_path_floyd(graph *g);
double graph_short_path_floyd_parallel(graph *g);
void printm(int **matrix, int cnt);

#endif
