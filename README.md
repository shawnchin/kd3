# Introduction

This is an implementation of a balanced 3D k-d tree used to speed up the selection of points within an interaction radius.

It is not meant to be a generic k-d tree solution. The existing APIs and functionality were designed based on the requirements of an exisiting application.

* The number of dimensions is fixed.
* Data points are defined as separate double arrays of equal size
* The tree is rebuilt often (points move) and searched very frequently (locating neighbours of each point).

Some optimisation were implemented, e.g. memory pooling of tree nodes and recycling of memory used by the search tree and iterator objects, but there is no doubt still more room for improvement.

# Usage example

```````C
#include "kd3/kdtree.h"
#define SEARCH_RADIUS 10
#define SEARCH_RADIUS_SQUARED 100

void some_function(void) {
  kdtree *tree = NULL;
  kdtree_iterator *result = NULL;
  
  double distance_squared, distance;
  double *x, *y, *z; /* array of points */
  size_t i, j; /* loop indices */

  x = malloc(sizeof(double) * SIZE);
  y = malloc(sizeof(double) * SIZE);
  z = malloc(sizeof(double) * SIZE);
  /* initialise values ..... */

  /* for each iteration */
  while(continue) {
    
    /* build tree based on current positions */
    kdtree_build(x, y, z, SIZE, &tree); /* tree obj recycled */
      
    /* for each point */
    for (i = 0; i < SIZE; i++) {

      /* search for neighbours */
      kdtree_search(tree, &result, x[i], y[i], z[i], SEARCH_RADIUS); /* result obj recycled */
      /* kdtree_iterator_sort(result); // if you need results in order */

      /* loop through each neighbour */
      j = kdtee_iterator_get_next(result);
      while (j != KDTREE_END) {
                
        /* second filter to ignore points beyond absolute distance */
        distance_squared = (x[i]-x[j])*(x[i]-x[j]) + 
                           (y[i]-y[j])*(y[i]-y[j]) + 
                           (z[i]-z[j])*(z[i]-z[j]);
        if (( i!= j) && (distance_squared <= SEARCH_RADIUS_SQUARED)) {
          distance = sqrt(distance_squared);
          /* do stuff with point i and neighbour j ... */
        }
        
        j = kdtee_iterator_get_next(result); /* get next */
      }
    }
        
    /* move points around ... */
    
    /* perform some additional operation for points that fall within a specific space */
    kdtree_search_space(tree, &result, -10.0, 10.0,  /* -10 <= x <= 10 */
                                       -30.3, 0.0,   /* -30.3 <= y <= 0 */
                                       -DBL_MAX, DBL_MAX);  /* Any z value */
    j = kdtee_iterator_get_next(result);
    while (j != KDTREE_END) {
      do_something_with_point(j);
      j = kdtee_iterator_get_next(result); /* get next */
    }
    
  }

  /* clean at the end */
  kdtree_delete(&tree);
  kdtree_iterator_delete(&result);
  free(x); free(y); free(z);
}
```````

Note that the search space for `kdtree_search()` is cube-shaped rather spherical. The last argument specifies the perpendicular distance between the center point and each face of the cube.

We leave the final filtering of points (discarding points that are beyond the search radius) to users as this involves calculating the absolute distance between points. Most use cases require the distance value within the inner loop anyway so it makes more sense to leave the calculation within the user code.

For a more generic search, you can also use `kdtree_search_space()` which searches within the space (3D box) defined by specifying the min and max values for each dimension.
