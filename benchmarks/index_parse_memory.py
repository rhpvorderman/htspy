import sys

from htspy.bam import BamIndex
from htspy.bgzf import VirtualFileOffset

import memory_profiler
import gc
import time


@memory_profiler.profile
def main():
    bai = BamIndex.from_file(sys.argv[1])
    for i, contig_index in enumerate(bai.contig_indices):
        bin_len = len(contig_index.bin_index)
        bin_size = sys.getsizeof(contig_index.bin_index)
        lin_len = len(contig_index.linear_index)
        lin_size = sys.getsizeof(contig_index.linear_index)
        if i < 40 or bin_len > 20:
            print(f"{i}\t{bin_len}\t{bin_size}\t{lin_len}\t{lin_size}")


if __name__ == "__main__":
    main()