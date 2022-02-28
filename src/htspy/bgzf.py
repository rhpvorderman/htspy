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
import typing
import zlib
from typing import Iterator, Optional

try:
    from isal import isal_zlib
except ImportError:
    isal_zlib = None  # type: ignore

GZIP_MAGIC = b"\x1f\x8b"
GZIP_MAGIC_INT = int.from_bytes(GZIP_MAGIC, "little", signed=False)

BGZF_MAX_BLOCK_SIZE = 0x10000  # 64K, 65536. Same as bgzf.h
BGZF_BLOCK_SIZE = 0xff00  # 65280. Same as bgzf.h

# XFL not set
BGZF_BASE_HEADER = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"  # noqa: E501


class BGZFError(IOError):
    pass


class VirtualFileOffset(typing.NamedTuple):
    coffset: int
    uoffset: int

    @classmethod
    def from_bytes(cls, b: bytes):
        virtual_offset, = struct.unpack("<Q", b)
        coffset = virtual_offset >> 16
        uoffset = virtual_offset & 0xffff
        return cls(coffset, uoffset)

    def to_integer(self):
        virtual_offset = self.coffset << 16
        virtual_offset &= self.uoffset
        return virtual_offset


class BGZFReader:
    def __init__(self, filename: str):
        self._file = open(filename, 'rb')
        self._buffer = io.BytesIO()
        self._buffer_size = 0
        if isal_zlib:
            self._decompress = isal_zlib.decompress
            self._crc32 = isal_zlib.crc32
        else:
            self._decompress = zlib.decompress
            self._crc32 = zlib.crc32  # type: ignore
        self._last_block_eof = False

    def close(self):
        self._buffer.close()
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def read_block(self):
        """Read the BGZF file up to the next block boundary"""
        if self._buffer_size:
            # There is some data in the buffer. Read it all and reset the buffer.
            block = self._buffer.read()
            self._buffer = io.BytesIO()
            self._buffer_size = 0
            if block:
                return block

        header = self._file.read(18)
        if not header:
            if self._last_block_eof:
                return b""
            else:
                raise EOFError("Truncated BGZF file. No EOF block found.")
        if len(header) < 18:
            raise EOFError(f"Truncated bgzf block.")
        magic, method, flags, mtime, xfl, os, xlen, si1, si2, slen, bsize = \
            struct.unpack("<HBBIBBHBBHH", header)
        if magic != GZIP_MAGIC_INT:
            raise BGZFError(f"Incorrect gzip magic: {struct.pack('<H', magic)}")
        if method != 8:  # Deflate method
            raise BGZFError(f"Unsupported compression method: {method}")
        if not flags & 4:
            raise BGZFError(f"Gzip block should contain an extra field.")
        if xlen < 6:
            raise BGZFError(f"XLEN too small")
        if not (si1 == 66 and si2 == 67 and slen == 2):
            raise BGZFError(f"Invalid BSIZE fields")
        # Skip other xtra fields.
        self._file.read(xlen - 6)
        block_size = bsize - xlen - 19
        block_peek = self._file.peek(1)
        if block_peek[0] == 1:  # No compression.
            length, inverse_length = struct.unpack("<xHH", self._file.read(5))
            if length != ~inverse_length & 0xFFFF or length != block_size - 5:
                raise BGZFError(f"Corrupted uncompressed block.")
            decompressed_block = self._file.read(length)
        else:
            block = self._file.read(block_size)
            if len(block) < block_size:
                raise EOFError(f"Truncated block.")
            # Decompress block, use the 64K as initial buffer size to avoid
            # resizing of the buffer. (Max block size before compressing is
            # slightly less than 64K for BGZF blocks). 64K is allocated faster
            # than sizes that are not powers of 2.
            decompressed_block = self._decompress(block,
                                                  wbits=-zlib.MAX_WBITS,
                                                  bufsize=65536)
        trailer = self._file.read(8)
        if len(trailer) < 8:
            raise EOFError(f"Truncated block.")
        crc, isize = struct.unpack("<II", trailer)
        if crc != self._crc32(decompressed_block):
            raise BGZFError("Checksum fail of decompressed block")
        if isize != len(decompressed_block):
            raise BGZFError("Incorrect length of decompressed blocks.")
        # Empty block signifies EOF.
        if not decompressed_block:
            self._last_block_eof = True
            # Jump back and check if there is more data. Empty blocks in
            # middle of the file should not return b"" but skipped.
            return self.read_block()
        self._last_block_eof = False
        return decompressed_block

    def __iter__(self):
        while True:
            block = self.read_block()
            if not block:
                return
            yield block

    def read(self, size=-1) -> bytes:
        if size == -1:
            return self.readall()
        current_pos = self._buffer.tell()
        if current_pos == self._buffer_size:
            # End of current buffer reached, read new block
            block = self.read_block()
            if not block:
                return b""
            self._buffer = io.BytesIO(block)
            self._buffer_size = len(block)
            current_pos = 0

        while size > (self._buffer_size - current_pos):
            block = self.read_block()
            if not block:
                break
            self._buffer.seek(self._buffer_size)
            self._buffer.write(block)
            self._buffer_size += len(block)
        self._buffer.seek(current_pos)
        return self._buffer.read(size)

    def readall(self) -> bytes:
        current_pos = self._buffer.tell()
        self._buffer.seek(self._buffer_size)  # move to end of buffer
        while True:
            block = self.read_block()
            if not block:
                break
            self._buffer.write(block)
        self._buffer.seek(current_pos)
        return self._buffer.read()


def _zlib_compress(data, level: int = -1, wbits: int = zlib.MAX_WBITS) -> bytes:
    """zlib.compress but with a wbits parameter."""
    compressobj = zlib.compressobj(level, wbits=wbits)
    return compressobj.compress(data) + compressobj.flush()


class BGZFWriter:
    def __init__(self, filename: str, compresslevel: Optional[int] = None):
        self._file = open(filename, 'wb')
        self._buffer = io.BytesIO(bytearray(BGZF_MAX_BLOCK_SIZE))
        self._buffer_size = 0
        self._buffer.seek(0)
        if isal_zlib:
            self._compress = isal_zlib.compress
            self._crc32 = isal_zlib.crc32
        else:
            self._compress = _zlib_compress
            self._crc32 = zlib.crc32  # type: ignore
        default_compresslevel = 1
        self.compresslevel = (compresslevel if compresslevel is not None
                              else default_compresslevel)

    def close(self):
        self.flush()
        self.write_eof_block()
        self._buffer.close()
        self._file.close()

    def write_eof_block(self):
        self._file.write(b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00"
                         b"\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00"
                         b"\x00\x00\x00\x00")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def flush(self):
        data_view = self._buffer.getbuffer()[:self._buffer_size]
        self.write_block(data_view)
        self._buffer_size = 0
        self._buffer.seek(0)

    def write_block(self, data):
        """Write a block of data immediately to the BGZF file as a block."""
        data_length = len(data)
        if data_length > BGZF_BLOCK_SIZE:
            raise ValueError(f"Cannot write data larger than "
                             f"{BGZF_BLOCK_SIZE} to a BGZF block.")
        self._file.write(BGZF_BASE_HEADER)
        if self.compresslevel:
            compressed_block = self._compress(data, self.compresslevel,
                                              wbits=-zlib.MAX_WBITS)
            # Length of the compressed block + generic gzip header (10 bytes) +
            # XLEN field (2 bytes)
            bgzf_block_size_bytes = struct.pack("<H", len(compressed_block) + 25)
            self._file.write(bgzf_block_size_bytes)
            self._file.write(compressed_block)
        else:
            size_and_deflate_header = struct.pack(
                "<HBHH",
                data_length + 30,  # +5 for deflate header, + 25 for rest of block.
                # Deflate block header: first bit signifying last block;
                # second and third bit 0 and 0 means uncompressed block.
                1,
                data_length,  # LEN
                ~data_length & 0xFFFF,  # NLEN
            )
            self._file.write(size_and_deflate_header)
            self._file.write(data)
        trailer = struct.pack("<II",
                              self._crc32(data),
                              data_length)
        self._file.write(trailer)

    def write(self, data):
        data_length = len(data)
        new_size = self._buffer_size + data_length
        if new_size > BGZF_BLOCK_SIZE:
            data_buffer = io.BytesIO(data)
            while data_buffer.tell() != data_length:
                written = self._buffer.write(
                    data_buffer.read(BGZF_BLOCK_SIZE - self._buffer_size))
                self._buffer_size += written
                if self._buffer_size == BGZF_BLOCK_SIZE:
                    self.flush()
        else:
            self._buffer.write(data)
            self._buffer_size = new_size
        return data_length
