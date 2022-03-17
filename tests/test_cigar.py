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

from htspy.bam import BAM_CBACK, BAM_CDEL, BAM_CDIFF, BAM_CEQUAL, \
    BAM_CHARD_CLIP, BAM_CINS, BAM_CMATCH, BAM_CPAD, BAM_CREF_SKIP, \
    BAM_CSOFT_CLIP, Cigar

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
    (1 << 4) | BAM_CMATCH,
    (20 << 4) | BAM_CINS,
    (300 << 4) | BAM_CDEL,
    (4000 << 4) | BAM_CREF_SKIP,
    (50_000 << 4) | BAM_CSOFT_CLIP,
    (600_000 << 4) | BAM_CHARD_CLIP,
    (7_000_000 << 4) | BAM_CPAD,
    (80_000_000 << 4) | BAM_CEQUAL,
    (268_435_435 << 4) | BAM_CDIFF,
    (9 << 4) | BAM_CBACK,
]

CIGAR_TUPLES = [
    (BAM_CMATCH, 1),
    (BAM_CINS, 20),
    (BAM_CDEL, 300),
    (BAM_CREF_SKIP, 4000),
    (BAM_CSOFT_CLIP, 50_000),
    (BAM_CHARD_CLIP, 600_000),
    (BAM_CPAD, 7_000_000),
    (BAM_CEQUAL, 80_000_000),
    (BAM_CDIFF, 268_435_435),
    (BAM_CBACK, 9)
]


def test_bam_cigar___init__():
    bam_cigar = Cigar(CIGAR_STRING)
    assert bam_cigar.number_of_operations == len(CIGAR_TUPLES)


def test_bam_cigar__str__():
    assert str(Cigar(CIGAR_STRING)) == CIGAR_STRING


def test_bam_cigar_to_tuples():
    assert list(Cigar(CIGAR_STRING)) == CIGAR_TUPLES


def test_bam_cigar_buffer():
    bam_cigar = Cigar(CIGAR_STRING)
    view = memoryview(bam_cigar)
    assert view.tobytes() == bam_cigar.raw
    assert view.nbytes == len(bam_cigar.raw)
    assert view.ndim == 1
    assert view.format == "I"
    assert view.itemsize == 4
    assert view.tolist() == CIGAR_NUMBER_LIST
    assert view.obj is bam_cigar
    assert len(view) == len(CIGAR_NUMBER_LIST)
    assert view[0] == CIGAR_NUMBER_LIST[0]
    assert view[8] == CIGAR_NUMBER_LIST[8]


def test_bam_cigar_from_bytes():
    bam_cigar = Cigar(CIGAR_STRING)
    assert Cigar.from_bytes(bam_cigar.raw) == bam_cigar


def test_bam_cigar_from_bytes_error():
    bam_cigar = Cigar(CIGAR_STRING)
    with pytest.raises(ValueError) as error:
        Cigar.from_bytes(bam_cigar.raw + b"12")
    error.match('Size of b must be a multiple of 4')


def test_bam_cigar_from_buffer():
    bam_cigar = Cigar(CIGAR_STRING)
    cigar_array = array.array("I", CIGAR_NUMBER_LIST)
    assert Cigar.from_buffer(cigar_array) == bam_cigar
