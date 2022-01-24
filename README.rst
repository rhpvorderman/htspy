pybam
=====

Fast pythonic BAM reading according to the `HTS format specifications
<http://samtools.github.io/hts-specs/>`_.

Differences with pysam
======================

This project was started to address some of pysam's less attractive features.
Due to pysam's long legacy, it is not possible to retroactively redesign pysam
without a massive break in backwards compatibility. Therefore pybam was
created.

+ Speed: pybam reads faster than pysam. (writing a work in progress as of this
  moment.)
+ Architecture: pybam interacts with BAM files and records directly, pysam only
  interacts with BAM files and records trough htslib.
+ Interface: pybam strives to offer a limited API that clearly communicates the
  limitations of the BAM format. Pysam has a very (very) extensive API with
  tons of features.

Design principles
=====================
+ Reading the `samv1 specification
  <https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf>`_
  (the BAM format section) should be a good introduction to this API.
+ The API must represent the underlying format as closely as possible. Example:
  BAM records do not store a reference name. Therefore BAM records may not
  have a ``reference_name`` attribute. (A method is acceptable).
+ Use Pythonic naming such as ``reference_id`` but refer to the name of the
  underlying attribute in the property documentation.
+ Do not use Python's properties to disguise from the user when intensive
  calculation is taking place. Example: ``cigar`` is an array of 32-bit
  integers each of which encodes for an operation and a number. A conversion to
  tuples or string requires quite some examples. Therefore ``cigar_as_string``
  or ``cigar_as_tuples`` must be methods. ``mate_unmapped`` on the other hand
  may be a property attribute as it translates directly to a bit set in the
  BAM record. Even though it takes a (simple) calculation to present that bit
  as a Python boolean.
+ Keep the API small. Having a gazillion methods on objects does not help
  with ease of programming.

Writing a pythonic BAM parser
=============================
The BAM format is very easy to parse. A clear specification is written. The
BAM format is designed such that whole parts can be copied into C-structs
directly on little-endian machines. Since each Python object is internally
also a C-struct, it is quite doable to write a BAM parser that creates Python
objects directly. This makes it much easier to manipulate BAM files than
letting htslib create a C-struct and then interacting with that C-struct
trough htslib. Pybam uses the Python C-API to create Python objects directly
from BAM files.

The BGZF format is a format that follows the gzip specification. It can be
handled in Python code entirely by using the ``struct`` and ``zlib`` modules
to handle the headers/trailers and compression/decompression respectively.
To speed up compression/decompression the `python bindings for ISA-L
<https://github.com/pycompression/python-isal>`_ were used.

