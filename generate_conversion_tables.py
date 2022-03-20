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


BASE_CODES = "=ACMGRSVTWYHKDBN"


def nucleotide_to_number_table():
    table = [str(-1) for _ in range(256)]
    for i, nuc in enumerate(BASE_CODES):
        table[ord(nuc)] = i
    return table


def number_to_nucleotide_table():
    table = ["" for _ in range(256)]
    for i, nuc in enumerate(BASE_CODES):
        for j, nuc2 in enumerate(BASE_CODES):
            index = (i << 4) | j
            bases_literal = f"{nuc}{nuc2}"
            bases_hex = f"0x{bases_literal.encode('ascii').hex()}"
            table[index] = bases_hex
    return table


def make_256_table(variable_name, table, row_size = 16):
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


def main():
    with open("src/htspy/_conversions.h", "wt", encoding="utf-8") as out:
        out.write('#include "stdint.h"\n')
        out.write('#include "htslib/sam.h"\n')
        out.write('\n')
        out.write(make_256_table(
            "static const char bam_cigar_table[256]",
            character_to_bam_op_table())
        )
        out.write('\n')
        out.write(make_256_table(
            "static const char nucleotide_to_number[256]",
            nucleotide_to_number_table()
        ))
        out.write('\n')
        out.write(make_256_table(
            "static const uint16_t number_to_nucleotide_pair[256]",
            number_to_nucleotide_table(),
            row_size=8
        ))


if __name__ == "__main__":
    main()
