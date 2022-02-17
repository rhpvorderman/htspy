implementation_ideas.rst
========================

Not yet implemented ideas notebook.

CIGAR
-----
Currently the raw CIGAR is stored as a bytes object. How to present this object
to the user?

+ Implement a BamCigar object that:
    + Uses a reference to the raw bytes object. (No copying required!)
    + Can be indexed 0-based.
    + Can be converted to string. (``__str__``)
    + Init method -> creation from string.
    + from_bytes method.
    + raw -> new reference to the underlying object.
    + Supports the buffer protocol.
    + Can be sliced.
    + Is immutable.

Tags
----
This is a hard problem. The current idea is to create a BamTag object. This
object has access to tag, value_type and value properties. But how to best
implement this object?