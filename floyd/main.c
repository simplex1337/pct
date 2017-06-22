#include <stdio.h>
#include "graph.h"

int main(int argc, char const *argv[])
{
    struct graph *g;
    g = graph_create(5);
    graph_set_edge(g, 1, 2, 10);
    graph_set_edge(g, 1, 4, 30);
    graph_set_edge(g, 1, 5, 100);
    graph_set_edge(g, 2, 3, 50);
    graph_set_edge(g, 3, 5, 10);
    graph_set_edge(g, 3, 4, 20);
    graph_set_edge(g, 4, 5, 60);
    struct graph *path;
    path = graph_short_path_floyd(g);
    printf("w =%d \n", graph_get_edge(path, 1, 2));
    printf("w =%d \n", graph_get_edge(path, 1, 4));
    printf("w =%d \n", graph_get_edge(path, 1, 5));
    printf("w =%d \n", graph_get_edge(path, 2, 3));
    printf("w =%d \n", graph_get_edge(path, 3, 5));
    printf("w =%d \n", graph_get_edge(path, 3, 4));
    printf("w =%d \n", graph_get_edge(path, 5, 5));
    printf("w =%d \n", graph_get_edge(path, 1, 3));
    graph_free(g);
    graph_free(path);
    return 0;
}