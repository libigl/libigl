.bf - bone forests
==================

------------------------------------------------------------------------

A .bf file contains a "bone forest". Normally a skeleton for linear blend skinning is a "bone tree" with a single root. But this format may store multiple trees, hence a forest.

Each line contains data about a vertex (joint) of the bone forest:

    [weight index] [parent index] [x] [y] [z] [undocument optional data]

Indices begin with 0. The weight index is -1 if the bone does not have an associated weight. The parent index is -1 for root nodes. The x,y,z coordinates are offset vectors from this joint's parent's location (for roots, an offset from the origin).
