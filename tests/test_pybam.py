import pytest

from pybam.bam import bam_iterator


@pytest.fixture
def raw_bam_record() -> bytes:
    # Seq consists of 7 AA characters
    seq = b'\x11\x11\x11\x10'
    # Quality is 35 for all bases
    quals = b'#######'
    l_seq = 7


    return b""

