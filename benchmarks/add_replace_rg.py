from htspy.bam import BamReader, BamWriter
import sys


if __name__ == "__main__":
    new_rg = sys.argv[2]
    with BamReader(sys.argv[1]) as bam_reader:
        with BamWriter(sys.argv[3],
                       bam_reader.header, compresslevel=0) as bam_writer:
            for record in bam_reader:
                record.set_tag("RG", new_rg)
                bam_writer.write(record)
