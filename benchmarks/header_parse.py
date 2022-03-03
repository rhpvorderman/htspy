import sys

from htspy.bam import BamReader

if __name__ == "__main__":
    with BamReader(sys.argv[1]) as bam:
        pass

