from pybam._bamrecord import bam_iterator
import gzip

with open("tests/no_hdr_sq_1.bam", "rb") as test_file:
    # Try to find start position
    test_data = gzip.decompress(test_file.read())

bamrecords = bam_iterator(test_data[268:])
first = next(bamrecords)
print(first.read_name)
