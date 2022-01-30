import pytest

from pybam._bamrecord import *

def test_bam_parsing() -> bytes:
    reference_id = -1
    pos = -1
    mapq = 99
    bin = 1001  # TODO get realistic value


    # Seq consists of 7 AA characters
    seq = b'\x11\x11\x11\x10'
    # Quality is 35 for all bases
    quals = b'#######'
    l_seq = 7

    return b""

