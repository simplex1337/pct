#include "graph.h"

#define min(x, y) (((x) < (y)) ? (x) : (y))

struct graph *graph_create(int nvertices)
{
    struct graph *g;
    g = malloc(sizeof(*g));
    g->nvertices = nvertices;
    g->visited = malloc(sizeof(int) * (nvertices + 1));
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

void graph_clear(struct graph *g)
{
    int i, j;
    for (i = 0; i < g->nvertices; i++) {
        g->visited[i] = 0;
        for (j = 0; j < g->nvertices; j++) {
            g->m[i][j] = INT_MAX;
        }
    }
    for (i = 0; i < g->nvertices; i++)
        g->m[i][i] = 0;
}

void graph_free(struct graph *g)
{
    for (int i = 0; i < g->nvertices; i++)
        free(g->m[i]);
    free(g->m);
    free(g);
}

void graph_set_edge(struct graph *g, int i, int j, int w)
{
    g->m[i - 1][j - 1] = w;
    g->m[j - 1][i - 1] = w;
}

int graph_get_edge(struct graph *g, int i, int j)
{
    return g->m[i - 1][j - 1];
}

int graph_nvertices(struct graph *g)
{
    return g->nvertices;
}

struct graph *graph_short_path_floyd(struct graph *g)
{
    struct graph *path;
    path = graph_create(g->nvertices);
    for (int i = 0; i < path->nvertices; i++)
        for (int j = 0; j < path->nvertices; j++)
            path->m[i][j] = g->m[i][j];
    for (int i = 0; i < path->nvertices; i++)
        for (int j = 0; j < path->nvertices; j++)
            for (int k = 0; k < path->nvertices; k++) {
                if(path->m[j][k] > path->m[j][i] + path->m[i][k])
                  path->m[j][k] = path->m[j][i] + path->m[i][k];
            }
    return path;
}
