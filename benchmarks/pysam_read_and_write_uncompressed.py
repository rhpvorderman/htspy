import sys

import pysam

with pysam.AlignmentFile(sys.argv[1], "rb") as bam_reader:
    with pysam.AlignmentFile(sys.argv[2], "wb0",
                             template=bam_reader) as bam_writer:
        for record in bam_reader:
            bam_writer.write(record)

