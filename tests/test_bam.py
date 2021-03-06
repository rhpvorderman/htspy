import array
import struct

from htspy._bam import BAM_CDIFF, BAM_CIGAR_SHIFT, BAM_CMATCH, \
    BAM_FUNMAP, BamRecord, Cigar, bam_iterator

import pytest


@pytest.fixture(scope="function")
def empty_bam() -> BamRecord:
    reference_id = -1
    pos = -1
    next_reference_id = -1
    next_pos = -1
    mapq = 0xff
    bin = 0  # Incorrect but roll with this for now
    flag = BAM_FUNMAP
    read_name = b""
    l_read_name = 1
    l_seq = 0
    seq = b''
    quals = b''
    cigar = b""
    n_cigar_op = 0
    tlen = 0
    tags = b''
    bam_struct = struct.pack("<iiBBHHHIiii",
                             reference_id, pos, l_read_name, mapq, bin,
                             n_cigar_op, flag, l_seq, next_reference_id,
                             next_pos, tlen)
    bam_record_without_block_size = (bam_struct + read_name + b"\x00" +
                                     cigar + seq + quals + tags)
    block_size = len(bam_record_without_block_size)
    bam_record = struct.pack("<I", block_size) + bam_record_without_block_size
    parsed_record = next(bam_iterator(bam_record))
    return parsed_record


def test_bam_parsing():
    reference_id = 3
    pos = 10000
    next_reference_id = -1
    next_pos = -1
    mapq = 99
    bin = 1001  # TODO get realistic value
    flag = 0
    read_name = b"my_forward_read/1"
    l_read_name = len(read_name) + 1
    l_seq = 7
    # Seq: GATTACA
    seq = b'\x41\x88\x12\x10'
    # Quality is 35 for all bases
    quals = b'#######'
    # 4 bases match, 3 bases differ.
    cigar = array.array("I", (BAM_CMATCH | (4 << BAM_CIGAR_SHIFT),
                              BAM_CDIFF | (3 << BAM_CIGAR_SHIFT)))
    n_cigar_op = len(cigar)
    tlen = 7  # Every base mapped.
    tags = b"RGZMySample\x00"
    bam_struct = struct.pack("<iiBBHHHIiii",
                             reference_id, pos, l_read_name, mapq, bin,
                             n_cigar_op, flag, l_seq, next_reference_id,
                             next_pos, tlen)
    bam_record_without_block_size = (bam_struct + read_name + b"\x00" +
                                     cigar.tobytes() + seq + quals + tags)
    block_size = len(bam_record_without_block_size)
    bam_record = struct.pack("<I", block_size) + bam_record_without_block_size
    parsed_record = next(bam_iterator(bam_record))
    assert parsed_record._block_size == block_size
    assert parsed_record._refID == reference_id
    assert parsed_record._pos == pos
    assert parsed_record._l_read_name == l_read_name
    assert parsed_record._mapq == mapq
    assert parsed_record._bin == bin
    assert parsed_record._n_cigar_op == n_cigar_op
    assert parsed_record._flag == flag
    assert parsed_record._l_seq == l_seq
    assert parsed_record._next_refID == next_reference_id
    assert parsed_record._next_pos == next_pos
    assert parsed_record._tlen == tlen
    assert parsed_record._read_name == read_name
    assert parsed_record._cigar == Cigar.from_buffer(cigar)
    assert parsed_record._seq == seq
    assert parsed_record._qual == quals
    assert parsed_record._tags == tags
    assert parsed_record.cigar == Cigar("4M3X")
    assert parsed_record.get_sequence() == "GATTACA"


def test_set_sequence_no_qual(empty_bam):
    empty_bam.set_sequence("GATTACA")
    assert empty_bam.get_sequence() == "GATTACA"
    assert empty_bam._seq == b'\x41\x88\x12\x10'
    assert empty_bam.qualities == b"\xff\xff\xff\xff\xff\xff\xff"


def test_set_sequence_with_qual(empty_bam):
    old_block_size = empty_bam._block_size
    empty_bam.set_sequence("GATTACA", b"\x1f\x1f\x1f\x1f\x1f\x1f\x1f")
    assert empty_bam.get_sequence() == "GATTACA"
    assert empty_bam.qualities == b"\x1f\x1f\x1f\x1f\x1f\x1f\x1f"
    assert empty_bam._block_size == \
           old_block_size + len(empty_bam._seq) + len(empty_bam.qualities)


def test_set_sequence_qual_wrong_type(empty_bam):
    with pytest.raises(TypeError) as error:
        empty_bam.set_sequence("GATTACA", "HFFFHF")
    assert error.match("qualities must be of type bytes")


def test_set_sequence_qual_wrong_length(empty_bam):
    with pytest.raises(ValueError) as error:
        empty_bam.set_sequence("GATTACA", b"FFFHF")
    assert error.match("sequence and qualities must have the same length")


def test_set_sequence_seq_wrong_type(empty_bam):
    with pytest.raises(TypeError) as error:
        empty_bam.set_sequence(b"GATTACA")
    assert error.match("sequence must be of type str")


def test_wrong_iupac_character_first_in_pair(empty_bam):
    with pytest.raises(ValueError) as error:
        empty_bam.set_sequence("XA")
    error.match("Not a IUPAC character: X")


def test_wrong_iupac_character_second_in_pair(empty_bam):
    with pytest.raises(ValueError) as error:
        empty_bam.set_sequence("AX")
    error.match("Not a IUPAC character: X")
