import array
import struct

from htspy._bam import BAM_CDIFF, BAM_CIGAR_SHIFT, BAM_CMATCH, \
    BAM_FUNMAP, Cigar, bam_iterator


def test_bam_parsing():
    reference_id = 3
    pos = 10000
    next_reference_id = -1
    next_pos = -1
    mapq = 99
    bin = 1001  # TODO get realistic value
    flag = BAM_FUNMAP
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
