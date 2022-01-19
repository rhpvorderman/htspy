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

GZIP_MAGIC = b"\x1f\x8b"
GZIP_MAGIC_INT = int.from_bytes(GZIP_MAGIC, "little", signed=False)


class BGZFError(IOError):
    pass


def decompress_bgzf_blocks(file: io.BufferedReader) -> Iterator[bytes]:
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
        block_size = bsize + 1
        block = file.read(block_size)
        if len(block) < block_size:
            raise EOFError(f"Truncated block at: {block_pos}")
        trailer = file.read(8)
        if len(trailer) < 8:
            raise EOFError(f"Truncated block at: {block_pos}")
        crc, isize = struct.unpack("<II", trailer)
        # Decompress block, use the isize as initial buffer size to avoid
        # resizing of the buffer.
        decompressed_block = zlib.decompress(block,
                                             wbits=-zlib.MAX_WBITS,
                                             bufsize=isize)
        if crc != zlib.crc32(decompressed_block):
            raise BGZFError("Checksum fail of decompressed block")
        if isize != len(decompressed_block):
            raise BGZFError("Incorrect length of decompressed blocks.")
        if decompressed_block == b"":
            # EOF Block found, check if we are at the EOF or if there is
            # another block.
            if not file.peek(1):
                return
        yield decompressed_block
