from htspy.bam import BamReader
import sys


if __name__ == "__main__":
    with BamReader(sys.argv[1]) as bam_reader:
        for record in bam_reader:
            seq = record.get_sequence()
            pass

