import io


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


def make_256_table(variable_name, table, row_size = 16):
    out = io.StringIO()
    out.write(variable_name + ' = {\n    ')
    i = 0
    for i, literal in enumerate(table):
        if i % row_size == 0 and i != 0:
            out.write(f" // {(i // row_size - 1) * row_size}-{i - 1}\n    ")
        out.write(literal.rjust(2, " ") + ", ")
    out.write(f" // {(i // row_size) * row_size}-{i}\n")
    return out.getvalue()


def main():
    with open("src/htspy/_conversions.h", "wt", encoding="utf-8") as out:
        out.write('#include "stdint.h"\n')
        out.write('#include "htslib/sam.h"\n')
        out.write('\n')
        out.write(make_256_table(
            "static const char bam_cigar_table[256]",
            character_to_bam_op_table())
        )
        out.write("};\n")


if __name__ == "__main__":
    main()
