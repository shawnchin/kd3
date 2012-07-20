#include <assert.h>
#include <stdio.h>
#include "kd3/kdtree.h"

double *x, *y, *z;

inline static void set_point(size_t idx, double X, double Y, double Z) {
  x[idx] = X; y[idx] = Y; z[idx] = Z;
}

static void initialise_points(void) {
  x = malloc(sizeof(double) * 11);
  y = malloc(sizeof(double) * 11);
  z = malloc(sizeof(double) * 11);
  
  set_point(0,  0.5,  0.5,  0.5);
  set_point(1,  0.5,  0.5,  0.5);
  set_point(2,  0.5,  0.5,  0.5);
  
  set_point(3,  0.0,  0.0,  0.0);
  set_point(4,  1.0,  0.0,  0.0);
  set_point(5,  1.0,  1.0,  0.0);
  set_point(6,  0.0,  1.0,  0.0);
  
  set_point(7,  0.0,  0.0,  1.0);
  set_point(8,  1.0,  0.0,  1.0);
  set_point(9,  1.0,  1.0,  1.0);
  set_point(10, 0.0,  1.0,  1.0);
}

static int cmp(const void* v1, const void* v2) {
  const int a = *((const int*)v1);
  const int b = *((const int*)v2);
  return ((a > b) ? 1 : ((a < b) ? -1 : 0));
}

/* expect v[] to be pre-sorted */
static void validate(kdtree_iterator *iter, size_t count, const size_t v[]) {
  size_t i = 0, *content;
  assert(iter != NULL);
  assert(iter->size == count);
  
  content = malloc(sizeof(size_t) * count);
  for (i = 0; i < count; i++) content[i] = kdtree_iterator_get_next(iter);
  assert(kdtree_iterator_get_next(iter) == KDTREE_END);
  
  qsort(content, count, sizeof(size_t), cmp);
  for (i = 0; i < count; i++) assert(content[i] == v[i]);
  
  free(content);
}

int main(void) {
  kdtree *tree = NULL;
  kdtree_iterator *iter = NULL;
  
  initialise_points();
  kdtree_build(x, y, z, 11, &tree);
  
  /* match none */
  kdtree_search(tree, &iter, -10, 0, 0, 9.999);
  const size_t e0[] = { 0 }; /* dummy value */
  validate(iter, 0, e0);
  
  /* match one */
  kdtree_search(tree, &iter, 0, 0, 0, 0.499);
  const size_t e1[] = { 3 };
  validate(iter, 1, e1);
  
  /* match all. intersect borders */
  kdtree_search(tree, &iter, 0.5, 0.5, 0.5, 0.5);
  const size_t e2[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  validate(iter, 11, e2);
  /* match all. beyond borders */
  kdtree_search(tree, &iter, 0.5, 0.5, 0.5, 100.0);
  validate(iter, 11, e2);
  
  /* front slice */
  kdtree_search(tree, &iter, 0.5, 0.5, 0.0, 0.5);
  const size_t e3[] = { 0, 1, 2, 3, 4, 5, 6 };
  validate(iter, 7, e3);
  
  /* back slice */
  kdtree_search(tree, &iter, 0.5, 0.5, 1.0, 0.5);
  const size_t e4[] = { 0, 1, 2, 7, 8, 9, 10 };
  validate(iter, 7, e4);
  
  /* using generic box search to search exactly top slice */
  kdtree_search_space(tree, &iter, 0.0, 1.0, 0.5, 1.0, 0.0, 1.0);
  const size_t e5[] = { 0, 1, 2, 5, 6, 9, 10 };
  validate(iter, 7, e5);
  
  printf("\n ---- ALL TESTS PASSED ---- \n");
  /* clean up */
  kdtree_iterator_delete(&iter);
  kdtree_delete(&tree);
  free(x); free(y); free(z);
  return 0;
}
