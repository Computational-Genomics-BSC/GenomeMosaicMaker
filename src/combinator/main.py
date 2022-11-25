# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import os
import sys
import argparse
import pysam
import bisect
import hashlib

# Add variant_extractor to PYTHONPATH
VARIANT_EXTRACTOR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '..', '..', 'dependencies', 'variant-extractor', 'src'))
sys.path.insert(0, VARIANT_EXTRACTOR_DIR)

from variant_extractor import VariantExtractor  # noqa


def _get_reads(file_obj, chrom, start, end):
    for read in file_obj.fetch(chrom, start, end):
        read_end = read.reference_start + read.query_length
        pair_end = read.next_reference_start + read.query_length
        if read.is_duplicate or read.is_secondary or read.is_supplementary or \
                read.next_reference_name != chrom or \
                pair_end > end or read.next_reference_start < start or \
                read_end > end or read.reference_start < start:
            continue
        yield read


def _open_output_file(output_file, template_file):
    # Get write mode based on the output file extension
    if output_file.endswith('.bam'):
        write_mode = 'wb'
    elif output_file.endswith('.sam'):
        write_mode = 'wh'
    elif output_file.endswith('.cram'):
        write_mode = 'wc'
    else:
        raise Exception('Invalid output file extension')

    # Open the output file
    return pysam.AlignmentFile(output_file, write_mode, template=template_file)


def _read_vcf(vcf_file, padding):
    # Create the VariantExtractor object
    variant_extractor = VariantExtractor(vcf_file)
    variants_df = variant_extractor.to_dataframe()

    # Create the dictionary of input files and zones
    zones = dict()
    input_files = dict()
    for _, row in variants_df.iterrows():
        var_obj = row['variant_record_obj']
        files = var_obj.info['FILES']
        if type(files) == str:
            files = [files]
        # Get absolute paths
        files = [os.path.abspath(f) for f in files]
        for idx, file in enumerate(files):
            if file not in input_files:
                input_files[file] = (idx, pysam.AlignmentFile(file))

        # Add the zone
        start_chrom = row['start_chrom']
        end_chrom = row['end_chrom']
        if start_chrom not in zones:
            zones[start_chrom] = []
        if end_chrom not in zones:
            zones[end_chrom] = []
        var_type = row['type_inferred']
        if var_type == 'SNV':
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files))
        elif var_type == 'INS':
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + row['length'] + padding, files))
        elif var_type == 'TRN':
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files))
            bisect.insort(zones[end_chrom], (row['end'] - padding, row['end'] + padding, files))
        else:
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['end'] + padding, files))

    return input_files, zones


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='Input VCF file')
    parser.add_argument('--outputs', '-o', required=True, nargs='+', help='Output alignment files')
    parser.add_argument('--padding', '-p', type=int, default=1000, help='Padding around the variants')

    args = parser.parse_args()

    # Convert everything to absolute paths
    args.input = os.path.abspath(args.input)
    args.outputs = [os.path.abspath(output) for output in args.outputs]

    input_files, zones = _read_vcf(args.input, args.padding)

    # Check for overlapping zones
    for chrom, chrom_zones in zones.items():
        for i in range(1, len(chrom_zones)):
            if chrom_zones[i][0] < chrom_zones[i - 1][1]:
                # Ask user if they want to continue
                print(
                    f'WARNING: Overlapping zones found in chromosome {chrom}. This may cause unrealistic final coverage.')
                print(f'         Zone {i - 1}: {chrom_zones[i - 1]}')
                print(f'         Zone {i}: {chrom_zones[i]}')
                print('Continue? [y/n]')
                answer = input()
                if answer != 'y':
                    sys.exit(1)

    # Create the output files
    output_files = []
    for idx, output_file in enumerate(args.outputs):
        output_files.append(_open_output_file(output_file, input_files[list(input_files.keys())[0]][1]))

    pending_unmapped_reads = dict()
    # Go through the zones and write the reads
    for chrom, chrom_zones in zones.items():
        for zone in chrom_zones:
            start, end, files = zone
            for file in files:
                file_idx, file_obj = input_files[file]
                for read in _get_reads(file_obj, chrom, start, end):
                    if read.mate_is_unmapped or read.is_unmapped:
                        # Check if the read is already in the pending reads
                        pending_mate = pending_unmapped_reads.pop(read.query_name, None)
                        if pending_mate is not None:
                            pending_unmapped_reads[read.query_name] = (read, file)
                    # Hash read name and file name to get an unique read name
                    read.query_name = hashlib.md5((f'{read.query_name}_{file}').encode()).hexdigest()
                    output_files[file_idx].write(read)

    # Write the missing unmapped reads
    for query_name, (mate_read, file) in pending_unmapped_reads.items():
        file_idx, file_obj = input_files[file]
        for read in file_obj.fetch(mate_read.next_reference_name, mate_read.next_reference_start):
            if read.is_unmapped and not read.is_duplicate and not read.is_secondary \
                    and not read.is_supplementary and read.query_name == query_name:
                read.query_name = hashlib.md5((f'{read.query_name}_{file}').encode()).hexdigest()
                output_files[file_idx].write(read)
                break

    # Close the output files
    for output_file in output_files:
        output_file.close()
