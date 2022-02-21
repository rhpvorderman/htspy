from htspy.bam import BamReader
import sys
import cProfile


def main():
    with BamReader(sys.argv[1]) as bam_reader:
        for record in bam_reader:
            pass


if __name__ == "__main__":
    cProfile.run("main()")
