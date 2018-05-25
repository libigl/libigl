.tgf - control handle graphs
============================

------------------------------------------------------------------------

A .tgf file contains a graph of describing a set of control handles/structures: point controls, bones of a skeleton and cages made of "cage edges".

The first part of the file consists of lines regarding each vertex of the graph. Each line reads:

    [index] [x] [y] [z] [undocument optional data]

Indices begin with 1 and should proceed in order. Then there should be a line with a sole:

    #

The next section concerns the edges of the graph. Each line corresponds to an edge:

    [source index] [dest index] [is bone] [is pseudo-edge] [is cage edge] [undocument other data]

Bone edges trump pseudo and cage edges.
