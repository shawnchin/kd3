/*!
 * \file kdtree.h
 *
 * \code
 *      Author: Shawn Chin
 *      Date  : July 2012
 *      Copyright (c) 2012 STFC Rutherford Appleton Laboratory
 * \endcode
 */
#include <stdlib.h> /* size_t */
#include <stdint.h> /* SIZE_MAX */

/* initial size for internal memory of iterator */
#define KDTREE_ITERATOR_INITIAL_SIZE 50

/* ratio to grow memory when iterator is full */
#define KDTREE_ITERATOR_GROWTH_RATIO 2

/* control value to indicate the end of iteration */
#ifndef SIZE_MAX
  #define KDTREE_END ((size_t)-1)
#else
  #define KDTREE_END SIZE_MAX
#endif

struct data_point {
  double x;
  double y;
  double z;
  size_t idx; /* index of original data point */
};

struct tree_node {
  struct tree_node *left;
  struct tree_node *right;
  double split;
  size_t idx;
};

struct boundaries {
  double min;
  double max;
};

struct space {
  struct boundaries dim[3];
};


typedef struct {
  size_t count;
  size_t max_nodes;
  size_t next_node;
  struct data_point *points;
  struct tree_node *node_data;
  struct tree_node *root;
} kdtree;

typedef struct {
  size_t *data;
  size_t capacity;
  size_t size;
  size_t current;
} kdtree_iterator;

void kdtree_build(double *x, double *y, double *z, size_t count, kdtree **tree);
void kdtree_delete(kdtree **tree_ptr);
void kdtree_search(kdtree *tree, kdtree_iterator **iter_ptr,
                   double x, double y, double z, double apothem);
void kdtree_search_space(kdtree *tree, kdtree_iterator **iter_ptr,
                         double x_min, double x_max,
                         double y_min, double y_max,
                         double z_min, double z_max);
size_t kdtree_iterator_get_next(kdtree_iterator *iter);
void kdtree_iterator_rewind(kdtree_iterator *iter);
void kdtree_iterator_sort(kdtree_iterator *iter);
void kdtree_iterator_delete(kdtree_iterator **iter_ptr);
