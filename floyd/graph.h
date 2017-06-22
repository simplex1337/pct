#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

struct g_path {
    int pathlen;
    int edge;
    int *path;
};

struct graph {
    int nvertices;              /* Число вершин */
    int **m;                    /* Матрица n x n */
    int *visited;
};

struct graph *graph_create(int nvertices);
void graph_clear(struct graph *g);
void graph_free(struct graph *g);
void graph_set_edge(struct graph *g, int i, int j, int w);
int graph_get_edge(struct graph *g, int i, int j);
int graph_nvertices(struct graph *g);
struct graph *graph_short_path_floyd(struct graph *g);


#endif
