.dmat - dense matrices
======================

------------------------------------------------------------------------

A .dmat file contains a dense matrix in column major order. It can contain ASCII or binary data. Note that it is uncompressed so binary only reduces the file size by 50%. But writing and reading binary is usually faster. In MATLAB, binary is almost 100x faster.

ASCII
-----

The first line is a header containing:

    [#cols] [#rows]

Then the coefficients are printed in column-major order separated by spaces.

Binary
------

Binary files will also contain the ascii header, but it should read:

    0 0

Then there should be another header containing the size of the binary part:

    [#cols] [#rows]

Then coefficients are written in column-major order in Little-endian 8-byte double precision IEEE floating point format.

**Note:** Line endings must be `'\n'` aka `char(10)` aka line feeds.
