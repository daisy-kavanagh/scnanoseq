#!/usr/bin/env python3

""" Given a fastq file and a blaze output file, this will extract the barcode
    and umi and place them in the header of the fastq, as well as stripping
    them from teh read.
"""
#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from concurrent.futures import ProcessPoolExecutor


def parse_args():
    parser = argparse.ArgumentParser(description="Extract barcode and UMI from reads using BLAZE output.")
    parser.add_argument("-i", "--input_file", required=True, help="Input FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("-b", "--barcode_file", required=True, help="BLAZE output file")
    parser.add_argument("-o", "--output_file", required=True, help="Output prefix")
    parser.add_argument("-f", "--barcode_format", required=True, help="Barcode/UMI format (e.g., 10X_3v3)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    return parser.parse_args()


def find_seq_indices(barcode, sequence, qualities):
    index = sequence.find(barcode)
    if index < 0:
        sequence = str(Seq(sequence).reverse_complement())
        qualities = qualities[::-1]
        index = sequence.find(barcode)
    return index, sequence, qualities


def strip_read_10X(bc_index, seq, quals, bc_len, umi_len, polyT_len):
    return {
        "bc": seq[bc_index:bc_index + bc_len],
        "bc_qual": quals[bc_index:bc_index + bc_len],
        "umi": seq[bc_index + bc_len:bc_index + bc_len + umi_len],
        "umi_qual": quals[bc_index + bc_len:bc_index + bc_len + umi_len],
        "r2_read": seq[bc_index + bc_len + umi_len + polyT_len:],
        "r2_qual": quals[bc_index + bc_len + umi_len + polyT_len:]
    }


def strip_read_by_format(bc_format, bc_index, seq, quals):
    if bc_format in ("10X_3v3", "10X_3v4", "10X_5v3"):
        return strip_read_10X(bc_index, seq, quals, 16, 12, 10)
    elif bc_format == "10X_5v2":
        return strip_read_10X(bc_index, seq, quals, 16, 10, 10)
    else:
        return None


def process_single_record(args):
    read_id, barcode, orig_seq, orig_quals, bc_format = args
    bc_index, seq, quals = find_seq_indices(barcode, orig_seq, orig_quals)
    if bc_index < 0:
        return None

    read_info = strip_read_by_format(bc_format, bc_index, seq, quals)
    if not read_info:
        return None

    bc_output = "\t".join([
        read_id,
        read_info["bc"],
        read_info["bc_qual"],
        read_info["umi"],
        read_info["umi_qual"]
    ])

    fq_output = "\n".join([
        "@" + read_id,
        read_info["r2_read"],
        "+",
        read_info["r2_qual"]
    ])

    return bc_output, fq_output


def extract_barcode(input_file, barcode_file, output_prefix, bc_format, threads):
    fastq_open = gzip.open if input_file.endswith(".gz") else open
    with fastq_open(input_file, "rt") as fq, \
         open(barcode_file, "r") as bc_in, \
         open(f"{output_prefix}.putative_bc_umi.tsv", "w") as bc_out, \
         open(f"{output_prefix}.fastq", "w") as fq_out:

        bc_out.write("read_id\tbc\tbc_qual\tumi\tumi_qual\n")

        chunk_size = 100000
        buffer = []

        records = SeqIO.parse(fq, "fastq")
        for bc_line, record in zip(bc_in, records):
            bc_line = bc_line.strip()
            try:
                _, barcode, *_ = bc_line.split(",")
            except ValueError:
                continue  # skip malformed line

            orig_seq = str(record.seq)
            orig_quals = "".join([chr(q + 33) for q in record.letter_annotations["phred_quality"]])
            buffer.append((record.id, barcode, orig_seq, orig_quals, bc_format))

            if len(buffer) >= chunk_size:
                with ProcessPoolExecutor(max_workers=threads) as executor:
                    for result in executor.map(process_single_record, buffer):
                        if result:
                            bc_out.write(result[0] + "\n")
                            fq_out.write(result[1] + "\n")
                buffer = []

        # Final chunk
        if buffer:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                for result in executor.map(process_single_record, buffer):
                    if result:
                        bc_out.write(result[0] + "\n")
                        fq_out.write(result[1] + "\n")


def main():
    args = parse_args()
    extract_barcode(args.input_file, args.barcode_file, args.output_file, args.barcode_format, args.threads)


if __name__ == "__main__":
    main()

