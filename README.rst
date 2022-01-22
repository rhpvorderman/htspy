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

Design decisions
=====================
The BAM format is very easy to parse. A clear specification is written. The
BAM format is designed such that whole parts can be copied into C-structs
directly on little-endian machines. Since each Python object is internally
also a C-struct, it is quite doable to write a BAM parser that creates Python
objects directly. This makes it much easier to manipulate BAM files than
letting htslib create a C-struct and then interacting with that C-struct
trough htslib.
