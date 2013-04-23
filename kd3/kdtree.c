/*!
 * \file kdtree.c
 *
 * \code
 *      Author: Shawn Chin
 *      Date  : July 2012
 *      Copyright (c) 2012 STFC Rutherford Appleton Laboratory
 * \endcode
 *
 * \brief Implementation of a balanced 3D k-d tree to speed up selection
 * of neighbouring points that are within a specific distance.
 *
 * This was not designed as a generic k-d tree solution. The APIs
 * and functionality was written based on the requirements of an
 * existing application.
 *  - Data points are defined in the application as separate double arrays
 *      of equal size (double *x, *y, *z)
 *  - The tree is rebuild frequently as the points move at the end of
 *      each iteration.
 *  - Searches are performed multiple times per iteration to locate points
 *      within a specific radius from the target
 *
 * Since the tree is to be rebuilt and search repeatedly, we reuse memory
 * where possible and avoid multiple allocations by provisioning memory
 * for tree nodes from a contiguous block of memory.
 *
 */
#include <assert.h> /* assert() */
#include <float.h>  /* DBL_MAX */
#include <string.h> /* memcpy() */
#include "kdtree.h"

/* dimensions hard coded to 3. Declare constants for convenience */
enum DIMENSIONS { DIM_X = 0, DIM_Y, DIM_Z, NDIMS };

/* Routines used for sorting points along different axes
 *
 * We've tried more complicate cmp routines which proceed
 * to compare the other axes in the event that two points
 * have the same value, however tests show that a basic
 * compare results in better overall performance.
 *
 */
#define CMP(v1,v2) ((v1 > v2) ? 1 : ((v1 < v2) ? -1 : 0))

static int cmp_x(const void *a1, const void *a2) {
  const struct data_point *A1 = (const struct data_point*)a1;
  const struct data_point *A2 = (const struct data_point*)a2;
  return CMP(A1->x, A2->x);
}

static int cmp_y(const void *a1, const void *a2) {
  const struct data_point *A1 = (const struct data_point*)a1;
  const struct data_point *A2 = (const struct data_point*)a2;
  return CMP(A1->y, A2->y);
}

static int cmp_z(const void *a1, const void *a2) {
  const struct data_point *A1 = (const struct data_point*)a1;
  const struct data_point *A2 = (const struct data_point*)a2;
  return CMP(A1->z, A2->z);
}

/* for sorting iterators */
static int cmp_size_t(const void *a1, const void *a2) {
  const size_t *A1 = (const size_t*)a1;
  const size_t *A2 = (const size_t*)a2;
  return CMP(*A1, *A2);
}

/* datatype for cmp function pointer */
typedef int(*cmp_func)(const void *, const void *);

/* static array of cmp functions indexed by dim */
const cmp_func func_select[] = { cmp_x, cmp_y, cmp_z };

/* declaration of internal functions */
inline static struct tree_node* _next_node(kdtree *tree);
inline static struct tree_node* _get_leaf_node(kdtree *tree, size_t offset);
inline static struct tree_node* _get_branch_node(kdtree *tree, double split);
static struct tree_node* _build_kdtree(size_t idx_from, size_t idx_to,
                                       size_t depth, kdtree *tree);
inline static kdtree_iterator* _iterator_new(void);
inline static void _iterator_reset(kdtree_iterator *iter);
inline static void _iterator_push(kdtree_iterator *iter, size_t value);
inline static void _explore_branch(kdtree *tree,
                                   struct tree_node *node,
                                   size_t depth,
                                   const struct space *search_space,
                                   const struct space *domain,
                                   kdtree_iterator *iter);
static void _search_kdtree(kdtree *tree,
                           struct tree_node *root,
                           size_t depth,
                           const struct space *search_space,
                           const struct space *domain,
                           kdtree_iterator *iter);
inline static int _point_in_search_space(const struct data_point *point,
                                         const struct space *search_space);
inline static int _completely_enclosed(const struct space *search_space,
                                       const struct space *domain);
inline static int _search_area_intersects(const struct space *search_space,
                                          const struct space *domain);
#ifndef _DEBUG_MODE
inline static int _is_leaf_node(const struct tree_node *node);
#else
#define _is_leaf_node(node) (node->left == NULL)
#endif

/* ---------------- Implementation of public APIs --------------------------- */

/* Build a 3D k-d tree based on the points stored in x, y, z arrays (with count
 * specifying the number of points).
 *
 * To optimise for cases where the data points may move and we need to rebuild
 * the tree for said points, we take in a reference to the kdtree object pointer
 * instead of simply returning the address of a new object.
 *
 * This allows the user to specify a NULL pointer when creating an new tree object,
 * or reuse the memory of the previously created object when rebuilding the tree
 * during the next iteration.
 *
 *   kdtree *tree = NULL;
 *   for (...) {
 *       kdtree_build(x, y, z, count, &tree);
 *   }
 *   kdtree_delete(&tree);
 *
 * Note that tree object can only be reused if the count is equal. Mismatching
 * counts will cause the previous object to be deleted and a new one built in
 * its place.
 *
 * To reduce the amount of checks, we do not handle cases where count < 0.
 * Do ensure that we're dealing with at lease two points
 *
 */
void kdtree_build(double *x, double *y, double *z, size_t count, kdtree **tree_ptr) {
  size_t i;
  kdtree *tree = *tree_ptr;
  
  /* sanity check */
  assert(count > 1);

  /* Reallocate object if ptr == NULL or if count does not match */
  if (!tree || tree->count != count) {
    if (tree) kdtree_delete(&tree); /* delete prev obj */
    
    /* allocate new object and update user's reference */
    tree = malloc(sizeof(kdtree));
    assert(tree != NULL);
    *tree_ptr = tree;
    
    /* initialise values and memory */
    tree->count = count;
    tree->max_nodes = ((count - 1) * 2) + 1;
    tree->points = malloc(sizeof(struct data_point) * count);
    tree->node_data = malloc(sizeof(struct tree_node) * tree->max_nodes);
    assert(tree->points != NULL);
    assert(tree->node_data != NULL);
  }

  /* reset control values */
  tree->next_node = 0;

  /* cache coordinates of each point and map to the idx of the point */
  for (i = 0; i < count; i++) {
    tree->points[i].idx = i;
    tree->points[i].x = x[i];
    tree->points[i].y = y[i];
    tree->points[i].z = z[i];
  }

  /* build tree and store ptr to root node */
  tree->root = _build_kdtree(0, count - 1, 0, tree);

}

/* search tree for points that fall within the 3d cube defined by
 * x, y, z, apothem where apothem is the distance from the point
 * to each side of the cube.
 */
void kdtree_search(kdtree *tree, kdtree_iterator **iter_ptr,
                   double x, double y, double z, double apothem) {
  assert(apothem >= 0.0);
  kdtree_search_space(tree, iter_ptr, 
                      x - apothem, x + apothem,
                      y - apothem, y + apothem,
                      z - apothem, z + apothem);
}

/* search tree for points that fall within the 3d box defined by
 * x_min, x_max, y_min, y_max, z_min, z_max.
 */
void kdtree_search_space(kdtree *tree, kdtree_iterator **iter_ptr,
                         double x_min, double x_max,
                         double y_min, double y_max,
                         double z_min, double z_max) {
  kdtree_iterator *iter = *iter_ptr;
  struct space search_space;
  struct space domain;

  /* sanity checks */
  assert(tree != NULL);
  
  /* The tree should have at least one point */
  assert(tree->root != NULL);
  assert(!_is_leaf_node(tree->root));

  /* Either create a new iterator or reset an exisiting one */
  if (iter != NULL) _iterator_reset(iter);
  else {
    iter = _iterator_new();
    *iter_ptr = iter; /* write back new ptr to obj */
  }

  /* define the search space */
  search_space.dim[DIM_X].min = x_min;
  search_space.dim[DIM_X].max = x_max;
  search_space.dim[DIM_Y].min = y_min;
  search_space.dim[DIM_Y].max = y_max;
  search_space.dim[DIM_Z].min = z_min;
  search_space.dim[DIM_Z].max = z_max;

  /* set initial domain to infinite space */
  domain.dim[DIM_X].min = -DBL_MAX;
  domain.dim[DIM_X].max =  DBL_MAX;
  domain.dim[DIM_Y].min = -DBL_MAX;
  domain.dim[DIM_Y].max =  DBL_MAX;
  domain.dim[DIM_Z].min = -DBL_MAX;
  domain.dim[DIM_Z].max =  DBL_MAX;

  /* search tree */
  _search_kdtree(tree, tree->root, 0, &search_space, &domain, iter);
}

/* Deallocates a tree object referenced by tree_ptr and sets the ptr to NULL */
void kdtree_delete(kdtree **tree_ptr) {
  kdtree *tree = *tree_ptr;
  if (tree == NULL) return;
  
  free(tree->points);
  free(tree->node_data);
  free(tree);
  *tree_ptr = NULL;
}

/* returns the next entry in the iteration, or KDTREE_END if the
 * end is reached */
size_t kdtree_iterator_get_next(kdtree_iterator *iter) {
  if (iter->current == iter->size) return KDTREE_END;
  return iter->data[iter->current++];
}

/* rewind the iterator */
void kdtree_iterator_rewind(kdtree_iterator *iter) {
  assert(iter != NULL);
  iter->current = 0;
}

/* deallocate memory associated with an iterator */
void kdtree_iterator_delete(kdtree_iterator **iter_ptr) {
  kdtree_iterator *iter = *iter_ptr;
  if (iter == NULL) return;
  
  free(iter->data);
  free(iter);
  *iter_ptr = NULL;
}

/* sort entries within the iterator */
void kdtree_iterator_sort(kdtree_iterator *iter) {
  qsort(iter->data, iter->size, sizeof(size_t), cmp_size_t);
}

/* --------------- INTERNAL ROUTINES ------------------------------- */


/* get pointer to the next available node within the node data cache */
inline static struct tree_node* _next_node(kdtree *tree) {
  assert(tree->next_node < tree->max_nodes);
  return &tree->node_data[tree->next_node++];
}

/* return a branch node */
inline static struct tree_node* _get_branch_node(kdtree *tree, double split) {
  struct tree_node *node = _next_node(tree);
  node->split = split;
  return node;
}

/* return a leaf node. Holds the index of the actual data point */
inline static struct tree_node* _get_leaf_node(kdtree *tree, size_t offset) {
  struct tree_node *node = _next_node(tree);
  node->left = NULL;
  node->right = NULL;
  node->idx = offset;
  return node;
}

#ifndef _DEBUG_MODE
/* determine if a node is a leaf node */
inline static int _is_leaf_node(const struct tree_node *node) {
  return ((node->left == NULL) && (node->right == NULL));
}
#endif

/* returns true if point is within search space */
inline static int _point_in_search_space(const struct data_point *point,
                                         const struct space *search_space) {
  return ((point->x <= search_space->dim[DIM_X].max) &&
          (point->x >= search_space->dim[DIM_X].min) &&

          (point->y <= search_space->dim[DIM_Y].max) &&
          (point->y >= search_space->dim[DIM_Y].min) &&
          
          (point->z <= search_space->dim[DIM_Z].max) &&
          (point->z >= search_space->dim[DIM_Z].min));
}

/* returns true if domain is completely enclosed within search space */
inline static int _completely_enclosed(const struct space *search_space,
                                       const struct space *domain) {
  return ((domain->dim[DIM_X].min <= search_space->dim[DIM_X].max) &&
          (domain->dim[DIM_X].min >= search_space->dim[DIM_X].min) &&
          (domain->dim[DIM_X].max <= search_space->dim[DIM_X].max) &&
          (domain->dim[DIM_X].max >= search_space->dim[DIM_X].min) &&
          
          (domain->dim[DIM_Y].min <= search_space->dim[DIM_Y].max) &&
          (domain->dim[DIM_Y].min >= search_space->dim[DIM_Y].min) &&
          (domain->dim[DIM_Y].max <= search_space->dim[DIM_Y].max) &&
          (domain->dim[DIM_Y].max >= search_space->dim[DIM_Y].min) &&
          
          (domain->dim[DIM_Z].min <= search_space->dim[DIM_Z].max) &&
          (domain->dim[DIM_Z].min >= search_space->dim[DIM_Z].min) &&
          (domain->dim[DIM_Z].max <= search_space->dim[DIM_Z].max) &&
          (domain->dim[DIM_Z].max >= search_space->dim[DIM_Z].min));
}


/* returns true if search space and domain
 *
 * It is easier to determine if two cubes are completely separate, so
 * we do just that and negate the return value.
 */
inline static int _search_area_intersects(const struct space *search_space,
                                          const struct space *domain) {
  return !((search_space->dim[DIM_X].min > domain->dim[DIM_X].max) ||
           (search_space->dim[DIM_X].max < domain->dim[DIM_X].min) ||
          
           (search_space->dim[DIM_Y].min > domain->dim[DIM_Y].max) ||
           (search_space->dim[DIM_Y].max < domain->dim[DIM_Y].min) ||
          
           (search_space->dim[DIM_Z].min > domain->dim[DIM_Z].max) ||
           (search_space->dim[DIM_Z].max < domain->dim[DIM_Z].min));
}

/* add all leaf nodes under a branch to the iterator */
static void _report_all_leaves(const kdtree *tree,
                               const struct tree_node *node,
                               kdtree_iterator *iter) {
  if (_is_leaf_node(node)) {
    _iterator_push(iter, tree->points[node->idx].idx);
  } else {
    _report_all_leaves(tree, node->left, iter);
    _report_all_leaves(tree, node->right, iter);
  }
}

/* convenience function to explore a sub-domain */
inline static void _explore_branch(kdtree *tree,
                                   struct tree_node *node,
                                   size_t depth,
                                   const struct space *search_space,
                                   const struct space *domain,
                                   kdtree_iterator *iter) {
  if (_is_leaf_node(node)) {
    if (_point_in_search_space(tree->points + node->idx, search_space)) {
    _iterator_push(iter, tree->points[node->idx].idx);
    }
  } else if (_search_area_intersects(search_space, domain)) {
    if (_completely_enclosed(search_space, domain)) {
      _report_all_leaves(tree, node, iter);
    } else {
      _search_kdtree(tree, node, depth + 1, search_space, domain, iter);
    }
  }
}

/* Recursively search the tree for points within a search space.
 * Results are appended to the iterator object.
 */
static void _search_kdtree(kdtree *tree,
                           struct tree_node *root,
                           size_t depth,
                           const struct space *search_space,
                           const struct space *domain,
                           kdtree_iterator *iter) {
  const size_t axis = depth % NDIMS;
  struct space new_domain;
  
  /* initialise boundaries for new domain */
  memcpy(&new_domain, domain, sizeof(struct space));
  
  /* explore left branch */
  new_domain.dim[axis].max = root->split;
  _explore_branch(tree, root->left, depth, search_space, &new_domain, iter);
  
  /* explore right branch */
  new_domain.dim[axis].max = domain->dim[axis].max; /* reset */
  new_domain.dim[axis].min = root->split;
  _explore_branch(tree, root->right, depth, search_space, &new_domain, iter);
}

/* internal routine to recursively build the kdtree */
static struct tree_node* _build_kdtree(size_t idx_from, size_t idx_to,
                                       size_t depth, kdtree *tree) {
  double split;
  struct tree_node *node;
  struct data_point *point;
  const size_t count = idx_to - idx_from + 1;
  const size_t mid   = idx_from + ((idx_to - idx_from) / 2);
  const size_t axis  = depth % NDIMS;

  /* if there is only one point, return a leaf node */
  if (count == 1) return _get_leaf_node(tree, idx_from);
  
  /* sort the points within this group to determine the median point
   * - This can be a potential performance bottleneck. There are methods
   *   to determine median in linear time, but that can get rather
   *   complicated. Will consider if this proves to be an issue.
  */
  qsort((void*)(tree->points + idx_from), count,
         sizeof(struct data_point), func_select[axis]);
  
  /* determine point where axis will be split */
  point = &tree->points[mid];
  split = (axis == 0) ? point->x : ((axis == 1) ? point->y : point->z);
  
  /* recursively build a tree for the left and right planes */
  node = _get_branch_node(tree, split);
  node->left  = _build_kdtree(idx_from, mid, depth + 1, tree);
  node->right = _build_kdtree(mid + 1, idx_to, depth + 1, tree);
  
  return node;
}

/* allocate and initialise a new iterator object */
inline static kdtree_iterator* _iterator_new(void) {
  kdtree_iterator *iter = malloc(sizeof(kdtree_iterator));
  assert(iter != NULL);
  
  iter->current = 0;
  iter->size = 0;
  iter->capacity = KDTREE_ITERATOR_INITIAL_SIZE;
  iter->data = malloc(sizeof(size_t) * iter->capacity);
  assert(iter->data != NULL);
  
  return iter;
}

/* resets and iterator so its memory can be reused */
inline static void _iterator_reset(kdtree_iterator *iter) {
  assert(iter != NULL);
  iter->size = 0;
  iter->current = 0;
}

/* add a new value into the iterator. Resize memory if full */
inline static void _iterator_push(kdtree_iterator *iter, size_t value) {
  if (iter->size == iter->capacity) { /* full. need to grow capacity */
    assert(KDTREE_ITERATOR_GROWTH_RATIO > 1.0);
    iter->capacity *= KDTREE_ITERATOR_GROWTH_RATIO;
    iter->data = realloc(iter->data, sizeof(size_t) * iter->capacity);
  }
  iter->data[iter->size++] = value;
}
