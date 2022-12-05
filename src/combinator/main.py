# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import os
import sys
import argparse
import pysam
import bisect
import hashlib
import threading, queue
from concurrent.futures import ThreadPoolExecutor

# Add variant_extractor to PYTHONPATH
VARIANT_EXTRACTOR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '..', '..', 'dependencies', 'variant-extractor', 'src'))
sys.path.insert(0, VARIANT_EXTRACTOR_DIR)

from variant_extractor import VariantExtractor  # noqa

RG_NAME = 'COMBINED'


def writer():
    # Call to_write.get() until it returns None
    for f, read in iter(to_write.get, None):
        f.write(read)

to_write = queue.Queue(maxsize=100000)
threading.Thread(target=writer).start()

def _is_valid_read(read):
    return not read.is_secondary and not read.is_supplementary


def _find_mate(read, original_query_name, file_obj):
    for mate_read in file_obj.fetch(read.next_reference_name, read.next_reference_start, read.next_reference_start + 1):
        if mate_read.query_name == original_query_name \
                and mate_read.is_read1 != read.is_read1 \
                and _is_valid_read(mate_read):
            return mate_read
    raise Exception(f'Mate not found for read {original_query_name} ({read.query_name}) in file {file_obj.filename}')

def _write_mates_from_file_read(mates_info):
    for info in mates_info:
        read, original_query_name, input_file_obj, output_file_obj = info
        mate = _find_mate(read, original_query_name, input_file_obj)
        mate.query_name = read.query_name
        mate.set_tag('RG', RG_NAME)
        to_write.put((output_file_obj, mate))

def _flush_pending_mates(pending_mates):
    # Group pending mates by input file
    mates_by_file = dict()
    for mate in pending_mates.values():
        file = mate[2]
        if file not in mates_by_file:
            mates_by_file[file] = []
        mates_by_file[file].append(mate)
    pending_mates.clear()

    pool = ThreadPoolExecutor()
    for _ in pool.map(_write_mates_from_file_read, mates_by_file.values()):
        pass
    pool.shutdown()


def _open_output_file(output_file, file_index, template_file):
    # Get write mode based on the output file extension
    if output_file.endswith('.bam'):
        write_mode = 'wb'
    elif output_file.endswith('.sam'):
        write_mode = 'wh'
    elif output_file.endswith('.cram'):
        write_mode = 'wc'
    else:
        raise Exception('Invalid output file extension')

    # Get the header from the template file
    header = template_file.header
    new_header = header.to_dict()
    # Remove the PG and RG lines
    new_header['PG'] = []
    new_header['RG'] = [{'ID': RG_NAME, 'SM': str(file_index)}]

    # Open the output file
    return pysam.AlignmentFile(output_file, write_mode, header=new_header)


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
                # input_files[file] = (idx, None)

        # Add the zone
        start_chrom = row['start_chrom']
        end_chrom = row['end_chrom']
        if start_chrom not in zones:
            zones[start_chrom] = []
        if end_chrom not in zones:
            zones[end_chrom] = []
        var_type = row['type_inferred']
        if var_type == 'SNV':
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files, str(var_obj)))
        elif var_type == 'INS':
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] +
                          row['length'] + padding, files, str(var_obj)))
        elif var_type == 'TRN':
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files, str(var_obj)))
            bisect.insort(zones[end_chrom], (row['end'] - padding, row['end'] + padding, files, str(var_obj)))
        elif var_type == 'INV':
            if row['start'] + padding > row['end'] - padding:
                bisect.insort(zones[start_chrom], (row['start'] - padding, row['end'] + padding, files, str(var_obj)))
            else:
                bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files, str(var_obj)))
                bisect.insort(zones[end_chrom], (row['end'] - padding, row['end'] + padding, files, str(var_obj)))
        else:
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['end'] + padding, files, str(var_obj)))
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

    print(f'Loaded {len(input_files)} input files')

    end_analysis = False
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
                else:
                    end_analysis = True
                    break
        if end_analysis:
            break

    print('Successfully checked for overlapping zones')
    # Create the output files
    output_files = []
    for idx, output_file in enumerate(args.outputs):
        output_files.append(_open_output_file(output_file, idx, input_files[list(input_files.keys())[0]][1]))

    pending_mates = dict()
    # Go through the zones and write the reads
    for chrom, chrom_zones in zones.items():
        for zone in chrom_zones:
            start, end, files = zone[:3]
            for filename in files:
                file_idx, input_file_obj = input_files[filename]
                output_file_obj = output_files[file_idx]
                for read in input_file_obj.fetch(chrom, start, end):
                    if not _is_valid_read(read):
                        continue
                    new_query_name = hashlib.md5((f'{read.query_name}_{filename}').encode()).hexdigest()
                    # Check if the read is already in the pending reads
                    pending_mate = pending_mates.pop(new_query_name, None)
                    if pending_mate is None:
                        pending_mates[new_query_name] = (read, read.query_name, input_file_obj, output_file_obj)
                    # Hash read name and file name to get an unique read name
                    read.query_name = new_query_name
                    read.set_tag('RG', RG_NAME)
                    to_write.put((output_file_obj, read))

    print(f'Flushing {len(pending_mates)} pending reads')
    # Write the missing mates
    _flush_pending_mates(pending_mates)

    # enqueue None to instruct the writer thread to exit
    to_write.put(None)

    # Close the output files
    for output_file in output_files:
        output_file.close()

    print('Done')