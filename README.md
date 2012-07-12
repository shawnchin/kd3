# Introduction

This is an implementation of a balanced 3D k-d tree used to speed up the selection of points within an interaction radius.

It is not meant to be a generic k-d tree solution. The existing APIs and functionality were designed based on the requirements of an exisiting application.

* The number of dimensions is fixed.
* Data points are defined as separate double arrays of equal size
* The tree is rebuilt often (points move) and searched very frequently (locating neighbours of each point).

Some optimisation were implemented, e.g. memory pooling of tree nodes and recycling of memory used by the search tree and iterator objects, but there is no doubt still more room for improvement.

# Usage example

<script src="https://gist.github.com/3097756.js"> </script>

Note that the search space is cube-shaped rather spherical. The last argument in `kdtree_search()` specifies the perpendicular distance between the center point and each face of the cube.

We leave the final filtering of points (discarding points that are beyond the search radius) to users as this involves calculating the absolute distance between points. Most use cases require the distance value within the inner loop anyway so it makes more sense to leave the calculation within the user code.
