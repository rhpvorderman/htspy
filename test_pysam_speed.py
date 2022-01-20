import pysam
import sys

if __name__ == "__main__":
    bam = pysam.AlignmentFile(sys.argv[1])
    for record in bam:
        pass
