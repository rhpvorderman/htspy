import cProfile
import sys

from htspy.bam import BamIndex


def main():
    BamIndex.from_file(sys.argv[1])


if __name__ == "__main__":
    cProfile.run("main()")