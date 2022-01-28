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
from typing import Iterator, List

from .bgzf import BGZFReader, BGZFWriter
from ._bamrecord import BamRecord, bam_iterator


class BAMFormatError(Exception):
    pass


class BamReader:
    def __init__(self, filename: str):
        self._file = BGZFReader(filename)
        self.header = b""
        self.references = []
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
        self.header = self._file.read(header_size)
        number_of_references, = struct.unpack("<I", self._file.read(4))
        for i in range(number_of_references):
            name_length, = struct.unpack("<I", self._file.read(4))
            name = self._file.read(name_length)
            seq_len, = struct.unpack("<I", self._file.read(4))
            self.references.append((name, seq_len))

    def __iter__(self) -> Iterator[BamRecord]:
        yield from bam_iterator(self._file.read_until_next_block())
        for block in iter(self._file):
            yield from bam_iterator(block)


class BamWriter:
    def __init__(self, filename: str, header: bytes, references: List[bytes]):
        self._file = BGZFWriter(filename)
        self._write_header(header, references)

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _write_header(self, header, references):
        self._file.write(b"BAM\x01")
        self._file.write(struct.pack("<I", len(header)))
        self._file.write(header)
        self._file.write(struct.pack("<I", len(references)))
        for name, seq_len in references:
            self._file.write(struct.pack("<I", len(name)))
            self._file.write(name)
            self._file.write(struct.pack("<I", seq_len))
        self._file.flush()

    def write(self, bam_record: BamRecord):
        self._file.write_block(bam_record.as_bytes())
