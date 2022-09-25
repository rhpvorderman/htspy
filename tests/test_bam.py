import array
import struct
import sys

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


def float32(f: float):
    """Represent float f with 32-bit precision."""
    reduced, = struct.unpack("<f", struct.pack("<f", f))
    return reduced


TEST_TAGS = (
        ("AB", b"ABAZ", "Z"),
        ("CD", b"CDZmystring\x00", "mystring"),
        ("EF", b"EFc" + struct.pack("<b", -20), -20),
        ("GH", b"GHC" + struct.pack("<B", 170), 170),
        ("IJ", b"IJs" + struct.pack("<h", -1024), -1024),
        ("KL", b"KLS" + struct.pack("<H", 65000), 65000),
        ("MN", b"MNi" + struct.pack("<i", -80_000), -80_000),
        ("OP", b"OPI" + struct.pack("<I", 3_000_000_000), 3_000_000_000),
        ("QR", b"QRf" + struct.pack("<f", 2.4),
            struct.unpack("<f", struct.pack("<f", 2.4))[0]),
        # B tag: first type letter, then count (of type uint32_t) then values.
        ("ST", b"STBc" + struct.pack("<Ibbb", 3, -20, 10, -126), [-20, 10, -126]),
        ("UV", b"UVBC" + struct.pack("<IBBB", 3, 65, 129, 203), [65, 129, 203]),
        ("WX", b"WXBs" + struct.pack("<Ihhh", 3, -4000, 4000, 2000), [-4000, 4000, 2000]),
        ("YZ", b"YZBS" + struct.pack("<IHHH", 3, 4000, 40000, 65000), [4000, 40000, 65000]),
        ("AA", b"AABi" + struct.pack("<Iiii", 3, -80000, 80000, 2_000_000_000), [-80000, 80000, 2_000_000_000]),
        ("BB", b"BBBI" + struct.pack("<IIII", 3, 80000, 2_000_000_000, 4_000_000_000), [80000, 2_000_000_000, 4_000_000_000]),
        ("CC", b"CCBf" + struct.pack("<Ifff", 3, 1.1, 2.2, 3.3), [float32(1.1), float32(2.2), float32(3.3)]),
        ("DD", b"DDBd" + struct.pack("<Iddd", 3, 1.1, 2.2, 3.3), [1.1, 2.2, 3.3]),
)


@pytest.mark.parametrize(["tag", "raw_tag", "result"], TEST_TAGS)
def test_get_tag(empty_bam, tag, raw_tag, result):
    empty_bam.tags = raw_tag
    ref_before = sys.getrefcount(raw_tag)
    retrieved_tag = empty_bam.get_tag(tag)
    if isinstance(retrieved_tag, memoryview):
        retrieved_tag = retrieved_tag.tolist()
    assert retrieved_tag == result
    del(result)
    assert sys.getrefcount(raw_tag) == ref_before


def concatenate_tags():
    for i, (tag, raw_tag, result) in enumerate(TEST_TAGS):
        yield tag, TEST_TAGS[i - 1][1] + raw_tag, result


@pytest.mark.parametrize(["tag", "raw_tag", "result"], concatenate_tags())
def test_get_tag_correct_skip(empty_bam, tag, raw_tag, result):
    empty_bam.tags = raw_tag
    retrieved_tag = empty_bam.get_tag(tag)
    if isinstance(retrieved_tag, memoryview):
        retrieved_tag = retrieved_tag.tolist()
    assert retrieved_tag == result


def truncated_tags():
    for tag, raw_tag, result in TEST_TAGS:
        # Checks tag_ptr to pyobject code
        for i in range(1, len(raw_tag)):
            yield tag, raw_tag[:-i]
        # checks skip tag code
        for i in range(1, len(raw_tag)):
            yield "XX", raw_tag[:-i]


@pytest.mark.parametrize(["tag", "trunc_tag"], truncated_tags())
def test_trucated_tag_error(empty_bam, tag, trunc_tag):
    empty_bam.tags = trunc_tag
    with pytest.raises(ValueError) as error:
        empty_bam.get_tag(tag)
    error.match("truncated tag")


@pytest.mark.parametrize(["tag", "raw_tag", "value"], TEST_TAGS)
def test_set_tag(empty_bam, tag, raw_tag, value):
    value_type = raw_tag[2:3].decode("ASCII")
    if value_type == "B":
        value_type = raw_tag[2:4].decode("ASCII")
        value = raw_tag[8:]
    empty_bam.set_tag(tag, value, value_type)
    assert empty_bam._tags == raw_tag


ARRAY_VALUES = [(array.array(python_type, [1, 2, 3]), bam_type)
                for python_type, bam_type in
                zip("bBhHiIf", "cCsSiIf")]


@pytest.mark.parametrize("value", ARRAY_VALUES)
def test_set_tag_array_format(value: array.array):
    value = value[0]
    empty_bam = BamRecord()
    empty_bam.set_tag("XX", value)
    arr: memoryview = empty_bam.get_tag("XX")
    assert list(arr) == list(value)
    assert arr.format == value.typecode


@pytest.mark.parametrize("value", ARRAY_VALUES[2:])
def test_set_tag_itemsize_error(value):
    value, typecode = value
    b = value.tobytes()[:-1]
    empty_bam = BamRecord()
    with pytest.raises(ValueError) as error:
        empty_bam.set_tag("XX", b, "B" + typecode)
    error.match("uffer size not a multiple of")


def test_set_tag_non_ascii_error():
    bam = BamRecord()
    with pytest.raises(ValueError) as error:
        bam.set_tag("Xë", "somestring")
    error.match("ASCII")


def test_set_tag_too_long_tag():
    bam = BamRecord()
    with pytest.raises(ValueError) as error:
        bam.set_tag("AVA", "somestring")
    error.match("2")
    error.match("length")


def test_set_tag_value_type_non_ascii():
    bam = BamRecord()
    with pytest.raises(ValueError) as error:
        bam.set_tag("AV", "somestring", "Č")
    error.match("ASCII")


def test_set_tag_value_type_too_long():
    bam = BamRecord()
    with pytest.raises(ValueError) as error:
        bam.set_tag("AV", "somestring", "ZZ")
    error.match("length")


@pytest.mark.parametrize(
    "value_type",
    ["c", "C", "s", "S", "i", "I", "f", "Z", "A",
     "Bc", "BC", "Bs", "BS", "Bi", "BI", "Bf" ]
)
def test_set_tag_wrong_types(value_type):
    # Strings do not support the buffer interface and are therefore
    # illegal for B type tags.
    value = "a" if value_type.startswith("B") else b"a"
    with pytest.raises(TypeError):
        BamRecord().set_tag("XX", value, value_type)


@pytest.mark.parametrize(["value", "value_type"], [
    ("String", "Z"),
    (1, "i"),
    (3.4, "f"),
    (array.array("b", [1]), "Bc"),
    (array.array("B", [1]), "BC"),
    (array.array("h", [1]), "Bs"),
    (array.array("H", [1]), "BS"),
    (array.array("i", [1]), "Bi"),
    (array.array("I", [1]), "BI"),
    (array.array("f", [1]), "Bf"),
])
def test_set_tag_autodetect_from_type(value, value_type):
    bam = BamRecord()
    bam.set_tag("XX", value)
    value_type_in_tag = bam._tags[2:3].decode("ASCII")
    if value_type_in_tag == "B":
        value_type_in_tag += bam._tags[3:4].decode("ASCII")
    assert value_type_in_tag == value_type


def test_set_tag_correctly_replaces():
    bam=BamRecord()
    # Test several tag replacements to make sure all code is run.
    bam.set_tag("XX", 1, 'C')  # Empty tag object gets filled with first value
    assert bam.get_tag("XX") == 1
    bam.set_tag("XY", 2, 'C')  # New tag gets added to existing tag object
    assert bam.get_tag("XY") == 2
    bam.set_tag("XZ", 3, 'C')  # Add another tag to have begin, middle and end tags.
    bam.set_tag("XY", 4, 'C')  # Replace middle tag in raw tag bytes object
    assert bam.get_tag("XY") == 4
    bam.set_tag("XX", 5, 'C')  # Replace first tag
    assert bam.get_tag("XX") == 5
    bam.set_tag("XX", 6, 'C')  # Replace last tag
    assert bam.get_tag("XX") == 6
