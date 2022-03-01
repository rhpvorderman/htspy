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
import os
import struct
import typing
from typing import BinaryIO, Dict, Iterator, List, Optional, Tuple

from ._bam import BamBlockBuffer as _BamBlockBuffer
from ._bam import BamRecord, bam_iterator
from .bgzf import BGZFReader, BGZFWriter, BGZF_BLOCK_SIZE, VirtualFileOffset


class BAMFormatError(Exception):
    pass


class BamReference(typing.NamedTuple):
    name: str
    length: int

    def to_bytes(self) -> bytes:
        nlen = struct.pack("<I", len(self.name))
        seq_len = struct.pack("<I", self.length)
        return nlen + self.name.encode('ascii', 'strict') + seq_len


class ContigIndex(typing.NamedTuple):
    binning_indices: Dict[int,
                          List[Tuple[VirtualFileOffset, VirtualFileOffset]]]
    linear_indices: List[VirtualFileOffset]
    reference_begin: Optional[VirtualFileOffset] = None
    reference_end: Optional[VirtualFileOffset] = None
    number_of_mapped_reads: Optional[int] = None
    number_of_unmapped_reads: Optional[int] = None

    @classmethod
    def from_fileobj(cls, fileobj: BinaryIO):
        binning_indices: Dict[
            int, List[Tuple[VirtualFileOffset, VirtualFileOffset]]] = {}
        linear_indices = []
        n_bin, = struct.unpack("<I", fileobj.read(4))
        for i in range(n_bin):
            bin, n_chunk = struct.unpack("<II", fileobj.read(8))
            chunk_list = []
            for j in range(n_chunk):
                chunk_list.append(
                    (VirtualFileOffset.from_bytes(fileobj.read(8)),
                     VirtualFileOffset.from_bytes(fileobj.read(8)))
                )
            binning_indices[bin] = chunk_list
        n_intv, = struct.unpack("<I", fileobj.read(4))
        for i in range(n_intv):
            linear_indices.append(
                VirtualFileOffset.from_bytes(fileobj.read(8)))
        pseudo_bins = binning_indices.get(37450)
        if pseudo_bins:
            assert len(pseudo_bins) == 2
            reference_begin, reference_end = pseudo_bins[0]
            nmapped_voffset, nunmapped_voffset = pseudo_bins[1]
            number_of_mapped_reads = nmapped_voffset.to_integer()
            number_of_unmapped_reads = nunmapped_voffset.to_integer()
            return cls(binning_indices, linear_indices, reference_begin,
                       reference_end, number_of_mapped_reads,
                       number_of_unmapped_reads)
        return cls(binning_indices, linear_indices)


class BamIndex(typing.NamedTuple):
    contig_indices: List[ContigIndex]
    number_of_unplaced_unmapped_reads: Optional[int] = None

    @classmethod
    def from_file(cls, filepath):
        with open(filepath, "rb") as file:
            magic = file.read(4)
            if magic != b"BAI\1":
                raise BAMFormatError(f"{filepath} is not a valid BAM index")
            n_ref, = struct.unpack("<I", file.read(4))
            contig_indices = [ContigIndex.from_fileobj(file)
                              for _ in range(n_ref)]
            n_no_coor_bytes = file.read(8)
            if n_no_coor_bytes:
                number_of_unplaced_unmapped_reads, = struct.unpack(
                    "<Q", n_no_coor_bytes)
                return cls(contig_indices, number_of_unplaced_unmapped_reads)
            return cls(contig_indices)


class BamHeader:
    def __init__(self, header: str, references=None):
        self.hd: Dict[str, str] = dict()
        self.sq: List[Dict[str, str]] = []
        self.rg: List[Dict[str, str]] = []
        self.pg: List[Dict[str, str]] = []
        self.co: List[str] = []
        self.references: List[BamReference] = []
        if references is not None:
            self.references = references[:]
        lines = header.splitlines(keepends=False)
        if lines[0].startswith("@HD\t"):
            _, self.hd = self.parse_tag_line(lines[0])
            self._check_tag_present("HD", "VN", self.hd)
            start_at = 1
        else:
            start_at = 0

        for line in lines[start_at:]:
            if line.startswith("@CO"):
                self.co.append(line.split("\t", 1)[1])
                continue
            record_type, tags_dict = self.parse_tag_line(line)
            if record_type == "SQ":
                self._check_tag_present("SQ", "SN", tags_dict)
                self._check_tag_present("SQ", "LN", tags_dict)
                self.sq.append(tags_dict)
            elif record_type == "RG":
                self._check_tag_present("RG", "ID", tags_dict)
                self.rg.append(tags_dict)
            elif record_type == "PG":
                self._check_tag_present("PG", "ID", tags_dict)
                self.pg.append(tags_dict)
            elif record_type == "HD":
                raise BAMFormatError(
                    "@HD must be the first line in the header")
            else:
                raise BAMFormatError(
                    f"Invalid record type in header: {record_type}")

    @staticmethod
    def parse_tag_line(line) -> Tuple[str, Dict[str, str]]:
        tags_dict = {}
        tags: List[str] = line.split("\t")
        record_type = tags.pop(0)
        record_type = record_type.lstrip("@")
        for tag in tags:
            # Tag can contain multiple colons. Only split on the first one.
            tag_name, tag_value = tag.split(":", maxsplit=1)
            tags_dict[tag_name] = tag_value
        return record_type, tags_dict

    @staticmethod
    def _check_tag_present(record_type, tag, tags_dict: Dict[str, str]):
        if tag not in tags_dict:
            raise BAMFormatError(
                f"{tag} is a mandatory tag on an @{record_type} line.")

    @staticmethod
    def _tag_dict_to_line(tag_dict):
        return "\t".join(f"{tag}:{value}" for tag, value in tag_dict.items())

    def to_sam_header(self) -> str:
        header_buffer = io.StringIO()
        if self.hd:  # Only write header line if not empty.
            header_buffer.write("@HD\t" + self._tag_dict_to_line(self.hd) + '\n')
        for tag_dict in self.sq:
            header_buffer.write("@SQ\t" + self._tag_dict_to_line(tag_dict) + '\n')
        for tag_dict in self.rg:
            header_buffer.write("@RG\t" + self._tag_dict_to_line(tag_dict) + '\n')
        for tag_dict in self.pg:
            header_buffer.write("@PG\t" + self._tag_dict_to_line(tag_dict) + '\n')
        for line in self.co:
            header_buffer.write("@CO\t" + line + '\n')
        return header_buffer.getvalue()

    def to_bytes(self) -> bytes:
        header_buffer = io.BytesIO()
        header_buffer.write(b"BAM\x01")
        sam_header = self.to_sam_header().encode('ascii')
        header_buffer.write(struct.pack("<I", len(sam_header)))
        header_buffer.write(sam_header)
        header_buffer.write(struct.pack("<I", len(self.references)))
        for reference in self.references:
            header_buffer.write(reference.to_bytes())
        return header_buffer.getvalue()


class BamReader:
    def __init__(self, filename: str, index_file: Optional[str] = None):
        self._file = BGZFReader(filename)
        if index_file is not None:
            self.index_file = index_file
        elif os.path.exists(filename + ".bai"):
            self.index_file = filename + ".bai"
        else:
            self.index_file = None
        self.index = None
        self.header: BamHeader
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
        sam_header = self._file.read(header_size)
        number_of_references, = struct.unpack("<I", self._file.read(4))
        references = []
        for i in range(number_of_references):
            name_length, = struct.unpack("<I", self._file.read(4))
            name = self._file.read(name_length)
            seq_len, = struct.unpack("<I", self._file.read(4))
            references.append(BamReference(name.decode('ascii'), seq_len))
        self.header = BamHeader(sam_header.decode('ascii'), references)

    def __iter__(self) -> Iterator[BamRecord]:
        for block in iter(self._file):
            yield from bam_iterator(block)


class BamWriter:
    def __init__(self, filename: str, header: BamHeader, compresslevel=None):
        self._file = BGZFWriter(filename, compresslevel)
        self.header = header
        self._write_header()
        self._buffer = _BamBlockBuffer(BGZF_BLOCK_SIZE)

    def close(self):
        self.flush()
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _write_header(self):
        self._file.write(self.header.to_bytes())
        self._file.flush()

    def flush(self):
        self._file.write_block(self._buffer.get_block_view())
        self._buffer.reset()

    def write(self, bam_record: BamRecord):
        if not self._buffer.write(bam_record):
            # Returned 0, buffer is full.
            self._file.write_block(self._buffer.get_block_view())
            self._buffer.reset()
            if not self._buffer.write(bam_record):
                # BamRecord to big for single block. Distribute over multiple
                # blocks.
                self._file.write(bam_record.to_bytes())
                self._file.flush()
