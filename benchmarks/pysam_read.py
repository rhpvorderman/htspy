import sys

import pysam

with pysam.AlignmentFile(sys.argv[1], "rb") as bam_reader:
    for record in bam_reader:
        pass
