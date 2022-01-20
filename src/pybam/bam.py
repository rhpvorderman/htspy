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

import struct
from typing import Iterator

from .bgzf import BGZFReader
from ._bamrecord import BamRecord, bam_iterator

class BAMFormatError(Exception):
    pass


class BamReader:
    def __init__(self, filename: str):
        self._file = BGZFReader(filename)
        self._read_header()

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _read_header(self):
        if self._file.read(4) != b"BAM\1":
            raise BAMFormatError("Not a BAM file")
        header_size, = struct.unpack("<I", self._file.read(4))
        self._file.read(header_size)  # Discard the header for now
        number_of_references, = struct.unpack("<I", self._file.read(4))
        for i in range(number_of_references):
            name_length, = struct.unpack("<I", self._file.read(4))
            self._file.read(name_length + 4)  # Discard name and size

    def __iter__(self):
        yield from bam_iterator(self._file.read_until_next_block())
        for block in iter(self._file):
            yield from bam_iterator(block)
