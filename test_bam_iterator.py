from pybam.bam import BamReader
import gzip

if __name__ == "__main__":
    with BamReader("tests/no_hdr_sq_1.bam") as bam_reader:
        for record in bam_reader:
            pass

