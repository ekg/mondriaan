#ifndef MATCHING_H_DONE

#define MATCHING_H_DONE

#include "solution.h"

#define TRUE 1
#define FALSE 0

void AugmentPath(int v, int *match, int *nmatch,
                 int *adj, int *start, int *length,
                 int *visited, int *vertex1, int *vertex2,
                 int *pred, int *dstart);

void add_vertex_to_graph(struct mat *a, struct solution *sol, int i,
                         unsigned char dir, unsigned char part);
void remove_vertex_from_graph(struct mat *a, struct solution *sol, int i,
                              unsigned char dir);
void update_graph_after_setting(struct mat *a, struct solution *sol, int i,
                                unsigned char dir, unsigned char part);
void update_graph_after_unsetting(struct mat *a, struct solution *sol, int i,
                                  unsigned char dir, unsigned char part);

#endif
