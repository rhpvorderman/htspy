import sys

from htspy.bam import BamIndex

if __name__ == "__main__":
    BamIndex.from_file(sys.argv[1])
