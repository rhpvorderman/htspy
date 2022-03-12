from htspy.bam import BamReader, BamWriter
import sys

import cProfile


def main():
    with BamReader(sys.argv[1]) as bam_reader:
        with BamWriter(sys.argv[2],
                       bam_reader.header, compresslevel=0) as bam_writer:
            for record in bam_reader:
                bam_writer.write(record)


if __name__ == "__main__":
    cProfile.run("main()")
