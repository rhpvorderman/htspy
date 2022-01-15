from pybam._bamrecord import bam_iterator
import zlib
with open("tests/no_hdr_sq_1.bam", "rb") as test_file:
    # Try to find start position
    test_data = zlib.decompress(test_file.read(), wbits=31)[280:]
print(test_data)

bamrecords = bam_iterator(test_data)
first = next(bamrecords)
