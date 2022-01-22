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

import io
import struct
import zlib
from typing import Iterator

try:
    from isal import isal_zlib
except ImportError:
    isal_zlib = None

GZIP_MAGIC = b"\x1f\x8b"
GZIP_MAGIC_INT = int.from_bytes(GZIP_MAGIC, "little", signed=False)


class BGZFError(IOError):
    pass


def decompress_bgzf_blocks(file: io.BufferedReader) -> Iterator[bytes]:
    if isal_zlib:
        decompress = isal_zlib.decompress
        crc32 = isal_zlib.crc32
    else:
        decompress = zlib.decompress
        crc32 = zlib.crc32
    while True:
        block_pos = file.tell()
        header = file.read(18)
        if len(header) < 18:
            raise EOFError(f"Truncated bgzf block at: {block_pos}")
        magic, method, flags, mtime, xfl, os, xlen, si1, si2, slen, bsize = \
            struct.unpack("<HBBIBBHBBHH", header)
        if magic != GZIP_MAGIC_INT:
            raise BGZFError(f"Invalid gzip block at: {block_pos}")
        if method != 8:  # Deflate method
            raise BGZFError(f"Unsupported compression method: {method} at"
                            f"block starting at: {block_pos}")
        if not flags & 4:
            raise BGZFError(f"Gzip block should contain an extra field. "
                            f"Block starts at: {block_pos}")
        if xlen < 6:
            raise BGZFError(f"XLEN too small at {block_pos}")
        if not (si1 == 66 and si2 == 67 and slen == 2):
            raise BGZFError(f"Invalid BSIZE fields at {block_pos}")
        # Skip other xtra fields.
        file.read(xlen - 6)
        block_size = bsize - xlen - 19
        block = file.read(block_size)
        if len(block) < block_size:
            raise EOFError(f"Truncated block at: {block_pos}")
        trailer = file.read(8)
        if len(trailer) < 8:
            raise EOFError(f"Truncated block at: {block_pos}")
        crc, isize = struct.unpack("<II", trailer)
        # Decompress block, use the isize as initial buffer size to avoid
        # resizing of the buffer.
        decompressed_block = decompress(block,
                                        wbits=-zlib.MAX_WBITS)
        if crc != crc32(decompressed_block):
            raise BGZFError("Checksum fail of decompressed block")
        if isize != len(decompressed_block):
            raise BGZFError("Incorrect length of decompressed blocks.")
        if decompressed_block == b"":
            # EOF Block found, check if we are at the EOF or if there is
            # another block.
            if not file.peek(1):
                return
        yield decompressed_block


class BGZFReader:
    def __init__(self, filename: str):
        self._file = open(filename, 'rb')
        self._block_iter = decompress_bgzf_blocks(self._file)  # type: ignore
        self._buffer = io.BytesIO()
        self._buffer_size = 0

    def close(self):
        self._buffer.close()
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self._block_iter

    def readall(self) -> bytes:
        current_pos = self._buffer.tell()
        self._buffer.seek(self._buffer_size) # move to end of buffer
        for block in self._block_iter:
            self._buffer.write(block)
        self._buffer.seek(current_pos)
        return self._buffer.read()

    def read(self, size=-1) -> bytes:
        if size == -1:
            return self.readall()
        current_pos = self._buffer.tell()
        if current_pos == self._buffer_size:
            # End of current buffer reached, read new block
            try:
                block = next(self._block_iter)
            except StopIteration:
                return b""
            self._buffer = io.BytesIO(block)
            self._buffer_size = len(block)
            current_pos = 0

        while size > (self._buffer_size - current_pos):
            try:
                block = next(self._block_iter)
            except StopIteration:
                break
            self._buffer.seek(self._buffer_size)
            self._buffer.write(block)
            self._buffer_size += len(block)
        self._buffer.seek(current_pos)
        return self._buffer.read(size)

    def read_until_next_block(self) -> bytes:
        """Read the BGZF file until the next BGZF block boundary."""
        if self._buffer.tell() == self._buffer_size:
            # Already at a block boundary, return next block
            try:
                return next(self._block_iter)
            except StopIteration:
                return b""
        # Otherwise, read the rest of the block in the buffer.
        return self._buffer.read()
