import sys

from htspy.bam import BamIndex

if __name__ == "__main__":
    bai = BamIndex.from_file(sys.argv[1])
    pass