import io
import struct


def character_to_bam_op_table():
    table = [str(-1) for _ in range(256)]
    table[ord("M")] = "BAM_CMATCH"
    table[ord("I")] = "BAM_CINS"
    table[ord("D")] = "BAM_CDEL"
    table[ord("N")] = "BAM_CREF_SKIP"
    table[ord("S")] = "BAM_CSOFT_CLIP"
    table[ord("H")] = "BAM_CHARD_CLIP"
    table[ord("P")] = "BAM_CPAD"
    table[ord("=")] = "BAM_CEQUAL"
    table[ord("X")] = "BAM_CDIFF"
    table[ord("B")] = "BAM_CBACK"
    return table


BASE_CODES = "=ACMGRSVTWYHKDBN"


def nucleotide_to_number_table():
    table = [str(-1) for _ in range(256)]
    for i, nuc in enumerate(BASE_CODES):
        table[ord(nuc)] = i
    # Make sure a NULL byte returns a NULL byte.
    table[0] = 0
    return table


def number_to_nucleotide_table_little_endian():
    table = ["" for _ in range(256)]
    for i, nuc in enumerate(BASE_CODES):
        for j, nuc2 in enumerate(BASE_CODES):
            index = (i << 4) | j
            bases_bytes = f"{nuc}{nuc2}".encode('ascii')
            base_integer, = struct.unpack("<H", bases_bytes)
            table[index] = hex(base_integer)
    return table


def make_table(variable_name, table, row_size = 16):
    out = io.StringIO()
    out.write(variable_name + ' = {\n    ')
    i = 0
    for i, literal in enumerate(table):
        if i % row_size == 0 and i != 0:
            out.write(f" // {(i // row_size - 1) * row_size}-{i - 1}\n    ")
        out.write(str(literal).rjust(2, " ") + ", ")
    out.write(f" // {(i // row_size) * row_size}-{i}\n")
    out.write("};\n")
    return out.getvalue()


def make_tag_value_function(sam_value_type,
                            c_type,
                            lower_bound,
                            upper_bound,
                            py_func,
                            tmp_type,
                            error_val):
    return (
        f"static int StorePyObjectValue_{sam_value_type}(PyObject *value, "
        f"void *value_store, uint8_t *tag)"
        f"{{\n"
        f"    {tmp_type} v = {py_func}(value);\n"
        f"    if ((v == {error_val}) && PyErr_Occurred()) {{\n"
        f"        return 0;\n"
        f"    }}\n"
        f"    if ((v < {lower_bound}) || (v > {upper_bound})) {{\n"
        f"        PyErr_Format(\n"
        f"            PyExc_ValueError,\n"
        f"            \"Tag '%c%c' with value_type '{sam_value_type}' should have a value \"\n"
        f"            \"between %ld and %ld.\",\n"
        f"            tag[0], tag[1], {lower_bound}, {upper_bound});\n"
        f"    }}\n"
        f"    (({c_type} *)value_store)[0] = ({c_type})v;\n"
        f"    return sizeof({c_type});\n"
        f"}}\n")


def main():
    with open("src/htspy/_conversions.h", "wt", encoding="utf-8") as out:
        out.write("#define PY_SSIZE_T_CLEAN\n")
        out.write("#include <Python.h>\n")
        out.write('#include "stdint.h"\n')
        out.write('#include "htslib/sam.h"\n')
        out.write('\n')
        out.write(make_table(
            "static const char bam_cigar_table[256]",
            character_to_bam_op_table())
        )
        out.write('\n')
        out.write(make_table(
            "static const char nucleotide_to_number[256]",
            nucleotide_to_number_table()
        ))
        out.write('\n')
        out.write(make_table(
            "static const uint16_t number_to_nucleotide_pair_le[256]",
            number_to_nucleotide_table_little_endian(),
            row_size=8
        ))


if __name__ == "__main__":
    main()
