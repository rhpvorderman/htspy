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

from typing import Iterable, Iterator, Tuple

class Cigar:
    def __init__(self, cigar_string: str): ...
    def __iter__(self) -> Iterator[Tuple[int, int]]: ...
    def __str__(self) -> str: ...
    def __eq__(self, other) -> bool: ...
    def __ne__(self, other) -> bool: ...

    @property
    def raw(self) -> bytes: ...
    @property
    def number_of_operations(self) -> int: ...

    @classmethod
    def from_bytes(cls, raw: bytes) -> Cigar: ...
    @classmethod
    def from_buffer(cls, buffer) -> Cigar: ...
    @classmethod
    def from_iter(cls, cigartuples: Iterable[Tuple[int, int]]) -> Cigar: ...


class BamRecord:
    _block_size: int
    _refID: int
    _pos: int
    _l_read_name: int
    _mapq: int
    _bin: int
    _n_cigar_op: int
    _flag: int
    _l_seq: int
    _next_refID: int
    _next_pos: int
    _tlen: int
    _read_name: bytes
    _cigar: bytes
    _seq: bytes
    _qual: bytes
    _tags: bytes

    reference_id: int
    position: int
    mapping_quality: int
    flag: int
    next_position: int
    template_length: int
    qualities: bytes

    @property
    def cigar(self) -> Cigar: ...

    @property
    def read_name(self) -> str: ...

    @property
    def tags(self) -> bytes: ...

    def to_bytes(self) -> bytes: ...


def bam_iterator(data) -> Iterator[BamRecord]: ...

class BamBlockBuffer:
    buffersize: int
    bytes_written: int

    def __init__(self, buffersize: int): ... 

    def write(self, bam_record: BamRecord) -> int: ... 

    def reset(self): ...  

    def get_block_view(self) -> memoryview: ...

BAM_CMATCH: int 
BAM_CINS: int 
BAM_CDEL: int 
BAM_CREF_SKIP: int 
BAM_CSOFT_CLIP: int 
BAM_CHARD_CLIP: int 
BAM_CPAD: int 
BAM_CEQUAL: int 
BAM_CDIFF: int 
BAM_CBACK: int
BAM_CIGAR_SHIFT: int

BAM_FPAIRED: int 
BAM_FPOPER_PAIR: int 
BAM_FUNMAP: int 
BAM_FREVERSE: int 
BAM_FMREVERSE: int 
BAM_FREAD1: int 
BAM_FREAD2: int 
BAM_FSECONDARY: int 
BAM_FQCFAIL: int 
BAM_FDUP: int 
BAM_FSUPPLEMENTARY: int 
