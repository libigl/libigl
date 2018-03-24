.rbr - ReAntTweakbar state file
===============================

------------------------------------------------------------------------

An .rbr file contains the saved values of the ReAntTweakBar class. It is used to load and save variables (and states specified via callbacks) stored in an AntTweakBar GUI.

Each line contains the name of the AntTweakBar item, the type of item and the value as a string:

    [name]: [type] [value]

As per AntTweakBar's own advice, names should not contain spaces. Names should also not contain colons (`:`). An example of a line looks like:

    my_rotation: TW_TYPE_QUAT4 0.0111272 -0.00101157 0.00648534 -0.999917

Not all AntTweakBar types are currently supported. See `igl/ReAntTweakbar.h` for an up-to-date list of supported types.
