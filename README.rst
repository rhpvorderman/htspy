htspy
=====

Fast pythonic BAM reading according to the `HTS format specifications
<http://samtools.github.io/hts-specs/>`_.

Just like `htsjdk <https://github.com/samtools/htsjdk>`_ htspy does not provide
bindings for htslib but natively implements the HTS specifications.

Goals:

+ Speed: due to implementing core features in C using the Python C-API htspy
  should be faster than any other python library supporting htslib formats
  (most notably `pysam <https://github.com/pysam-developers/pysam>`_).
+ Pythonic: features should be easy to use and feel intuitive for
  experienced Python programmers.
+ Architecture: interact natively with HTS formats rather than trough another
  library.

Currently the focus is on BAM, SAM and BGZF files. For FASTQ files there is
already an excellent library: `dnaio <https://www.github.com/marcelm/dnaio>`_.
For VCF there is already another very fast library: `cyvcf2
<https://https://github.com/brentp/cyvcf2>`_.

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
  tuples or string requires quite some calculations. Therefore ``cigar_as_string``
  or ``cigar_as_tuples`` must be methods. ``mate_unmapped`` on the other hand
  may be a property attribute as it translates directly to a bit set in the
  BAM record. It does not require calculation (merely a check if the bit is set).
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

