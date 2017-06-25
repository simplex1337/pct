#include "graph.h"
#include <string.h>

int main()
{
    graph *g;
    int n = 1000;
    g = graph_create(n);
    for (int i = 1; i < n; i++) {
        graph_set_edge(g, i, n, 4);
    }
    /*graph_set_edge(g, 1, 2, 10);
       graph_set_edge(g, 1, 4, 30);
       graph_set_edge(g, 1, 5, 100);
       graph_set_edge(g, 2, 3, 50);
       graph_set_edge(g, 3, 5, 10);
       graph_set_edge(g, 3, 4, 20);
       graph_set_edge(g, 4, 5, 60); */
    graph *g_path, *g_path1;
    g_path = graph_create(g->nvertices);
    g_path1 = graph_create(g->nvertices);
    for (int i = 0; i < g->nvertices; i++)
        for (int j = 0; j < g->nvertices; j++) {
            g_path->m[i][j] = g->m[i][j];
            g_path1->m[i][j] = g->m[i][j];
        }
    // printm(g_path->m, g_path->nvertices);
    double t = graph_short_path_floyd_serial(g_path);
    double t1 = graph_short_path_floyd_parallel(g_path1);
    // printm(g_path->m, g_path->nvertices);
    printf("w = %d\n", graph_get_edge(g_path, 1, n));
    printf("Time ser (sec): %.6f\nTime parallel (sec): %.6f\nS(n): %.6f\n", t, t1, t / t1);
    graph_free(g);
    graph_free(g_path);
    graph_free(g_path1);
    return 0;
}
