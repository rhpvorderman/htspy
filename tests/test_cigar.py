# Copyright (c) 2022 Ruben Vorderman
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import array

from htspy.bam import BamCigar, CigarOp

import pytest

CIGAR_STRING = ("1M"  # Test different digit lengths
                "20I"
                "300D"
                "4000N"
                "50000S"
                "600000H"
                "7000000P"
                "80000000="
                "268435435X"  # Maximum value: 2^28 - 1
                "9B")

CIGAR_NUMBER_LIST = [
    (1 << 4) | CigarOp.MATCH,
    (20 << 4) | CigarOp.INS,
    (300 << 4) | CigarOp.DEL,
    (4000 << 4) | CigarOp.REF_SKIP,
    (50_000 << 4) | CigarOp.SOFT_CLIP,
    (600_000 << 4) | CigarOp.HARD_CLIP,
    (7_000_000 << 4) | CigarOp.PAD,
    (80_000_000 << 4) | CigarOp.EQUAL,
    (268_435_435 << 4) | CigarOp.DIFF,
    (9 << 4) | CigarOp.BACK,
]

CIGAR_TUPLES = [
    (CigarOp.MATCH, 1),
    (CigarOp.INS, 20),
    (CigarOp.DEL, 300),
    (CigarOp.REF_SKIP, 4000),
    (CigarOp.SOFT_CLIP, 50_0000),
    (CigarOp.HARD_CLIP, 600_000),
    (CigarOp.PAD, 7_000_000),
    (CigarOp.EQUAL, 80_000_000),
    (CigarOp.DIFF, 268_435_435),
    (CigarOp.BACK, 9)
]

def test_bam_cigar___init__():
    bam_cigar = BamCigar(CIGAR_STRING)
    assert bam_cigar.number_of_operations == len(CIGAR_TUPLES)

def test_bam_cigar__str__(bam_cigar):
    assert str(BamCigar(CIGAR_STRING)) == CIGAR_STRING


def test_bam_cigar_to_tuples():
    assert list(BamCigar(CIGAR_STRING)) == CIGAR_TUPLES



def test_bam_cigar_buffer(bam_cigar):
    bam_cigar = BamCigar(CIGAR_STRING)
    view = memoryview(bam_cigar)
    assert view.tobytes() == bam_cigar.raw
    assert view.nbytes == len(bam_cigar.raw)
    assert view.ndim == 1
    assert view.format == "I"
    assert view.itemsize == 4
    assert view.tolist() == CIGAR_NUMBER_LIST
    assert view.obj is bam_cigar
    assert len(view) == 9
    assert view[0] == CIGAR_NUMBER_LIST[0]
    assert view[8] == CIGAR_NUMBER_LIST[8]


def test_bam_cigar_from_bytes(bam_cigar):
    assert BamCigar.from_bytes(bam_cigar.raw) == bam_cigar


def test_bam_cigar_from_bytes_error(bam_cigar):
    with pytest.raises(ValueError) as error:
        BamCigar.from_bytes(bam_cigar.raw + b"12")
    error.match("bytes length not a multiple of 4")


def test_bam_cigar_from_buffer(bam_cigar):
    cigar_array = array.array("I", CIGAR_NUMBER_LIST)
    assert BamCigar.from_buffer(cigar_array) == bam_cigar


def test_bam_cigar_from_buffer_incompatible_format(bam_cigar):
    cigar_array = array.array("i", CIGAR_NUMBER_LIST[:8])
    with pytest.raises(BufferError) as error:
        BamCigar.from_buffer(cigar_array)
    error.match("Invalid buffer format, expected: 'I', got: 'i'")
