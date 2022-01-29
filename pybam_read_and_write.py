from pybam.bam import BamReader, BamWriter
import sys


if __name__ == "__main__":
    with BamReader(sys.argv[1]) as bam_reader:
        with BamWriter(sys.argv[2],
                       bam_reader.header) as bam_writer:
            for record in bam_reader:
                bam_writer.write(record)
