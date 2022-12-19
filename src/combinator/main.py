# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import os
import sys
import argparse
import pysam
import bisect
import hashlib
from concurrent.futures import ProcessPoolExecutor
import logging

# Add variant_extractor to PYTHONPATH
VARIANT_EXTRACTOR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '..', '..', 'dependencies', 'variant-extractor', 'src'))
sys.path.insert(0, VARIANT_EXTRACTOR_DIR)

from variant_extractor import VariantExtractor  # noqa

RG_NAME = 'COMBINED'


def _write_read(output_file_obj, read, input_filename):
    new_query_name = hashlib.md5((f'{read.query_name}_{input_filename}').encode()).hexdigest()
    read.query_name = new_query_name
    read.set_tag('RG', RG_NAME)
    output_file_obj.write(read)


def _write_reads_in_zones(zones, input_filename, input_file_obj, output_file_obj):
    pending_mates = dict()
    for chrom, chrom_zones in zones.items():
        for zone in chrom_zones:
            for read in input_file_obj.fetch(chrom, zone[0], zone[1]):
                if read.is_secondary or read.is_supplementary:
                    continue
                mate_read = pending_mates.pop(read.query_name, None)
                if mate_read is None:
                    pending_mates[read.query_name] = (read, read.query_name)
                else:
                    _write_read(output_file_obj, read, input_filename)
                    _write_read(output_file_obj, mate_read[0], input_filename)
    return pending_mates


def _find_mate(pending_mates, already_found_set, read, original_query_name, file_obj):
    found_mates = []
    for mate_read in file_obj.fetch(read.next_reference_name, read.next_reference_start, read.next_reference_start + 1):
        if mate_read.is_secondary or mate_read.is_supplementary:
            continue
        pending_read, _ = pending_mates.get(mate_read.query_name, (None, None))
        if pending_read is None or mate_read.is_read1 == pending_read.is_read1:
            continue
        # Check if the mate is the one we are looking for
        if mate_read.query_name == original_query_name:
            found_mates.append(mate_read)
            return found_mates
        elif mate_read.query_name not in already_found_set:
            # This read is another pending mate
            found_mates.append(mate_read)
    logging.warning(f'Mate not found for read {original_query_name} ({read.query_name}) in file {file_obj.filename}. Ignoring it.')
    return found_mates


def _write_pending_mates(pending_mates, input_file_obj, output_file_obj, input_filename):
    already_found_set = set()
    for mate_read, original_query_name in pending_mates.values():
        if original_query_name in already_found_set:
            continue
        reads = _find_mate(pending_mates, already_found_set, mate_read, original_query_name, input_file_obj)
        for read in reads:
            already_found_set.add(read.query_name)
            original_read = pending_mates.get(read.query_name)[0]
            _write_read(output_file_obj, original_read, input_filename)
            _write_read(output_file_obj, read, input_filename)


def _write_zones(zones, file_index, input_filename, output_filename, multithreading, max_memory, fasta_ref):
    input_file_obj = pysam.AlignmentFile(
        input_filename, threads=2 if multithreading else 1, reference_filename=fasta_ref)
    output_filename_unsorted = f'{output_filename}.unsorted.{output_filename.split(".")[-1]}'
    output_file_obj = _open_output_file(output_filename_unsorted, file_index, input_file_obj, multithreading, fasta_ref)

    # Write reads
    pending_reads = _write_reads_in_zones(zones, input_filename, input_file_obj, output_file_obj)
    logging.debug(f'Found {len(pending_reads)} pending reads in {input_filename} > {output_filename}')

    # Write pending reads
    _write_pending_mates(pending_reads, input_file_obj, output_file_obj, input_filename)

    output_file_obj.close()
    input_file_obj.close()

    logging.debug(f'Finished reading {input_filename} > {output_filename}')

    max_memory_mb = int(max_memory * 1024)
    # Sort the output file
    pysam.sort('-o', output_filename, '-m', f'{max_memory_mb}M', '-T', f'{output_filename}_temp',
               output_filename_unsorted, '-@', '2' if multithreading else '1')
    os.remove(output_filename_unsorted)

    logging.debug(f'Finished sorting {output_filename}')


def _open_output_file(output_file, file_index, template_file, multithreading, fasta_ref):
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
    return pysam.AlignmentFile(output_file, write_mode, header=new_header, threads=2 if multithreading else 1, reference_filename=fasta_ref)


def _merge_files(processes, output_files, output_file):
    pysam.merge('-f', '-@', str(processes), *output_files, '-o', output_file)
    for output_file in output_files:
        os.remove(output_file)


def _read_vcf(vcf_file, padding):
    # Create the VariantExtractor object
    variant_extractor = VariantExtractor(vcf_file)
    variants_df = variant_extractor.to_dataframe()

    # Create the dictionary of input files and zones
    zones = dict()
    for _, row in variants_df.iterrows():
        var_obj = row['variant_record_obj']
        files = var_obj.info['FILES']
        if type(files) == str:
            files = [files]
        # Get absolute paths
        files = [os.path.abspath(f) for f in files]
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
    return zones


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='Input VCF file')
    parser.add_argument('--outputs', '-o', required=True, nargs='+', help='Output alignment files')
    parser.add_argument('--padding', '-p', type=int, default=1000, help='Padding around the variants')
    parser.add_argument('--maximum-processes', '-mp', type=int, default=1, help='Maximum number of processes to use')
    parser.add_argument('--maximum-memory', '-mm', type=float, default=32, help='Maximum memory to use (in GiB)')
    parser.add_argument('--multithreading', '-t', action='store_true', help='Use multithreading')
    parser.add_argument('--fasta-ref', '-f', type=str, help='Fasta reference file (used for CRAM files)')

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG)

    # Convert everything to absolute paths
    args.input = os.path.abspath(args.input)
    args.outputs = [os.path.abspath(output) for output in args.outputs]

    zones = _read_vcf(args.input, args.padding)

    logging.debug(f'Loaded zones')

    end_analysis = False
    # Check for overlapping zones
    for chrom, chrom_zones in zones.items():
        for i in range(1, len(chrom_zones)):
            if chrom_zones[i][0] < chrom_zones[i - 1][1]:
                # Ask user if they want to continue
                logging.debug(
                    f'WARNING: Overlapping zones found in chromosome {chrom}. This may cause unrealistic final coverage.')
                logging.debug(f'         Zone {i - 1}: {chrom_zones[i - 1]}')
                logging.debug(f'         Zone {i}: {chrom_zones[i]}')
                logging.debug('Continue? [y/n]')
                answer = input()
                if answer != 'y':
                    sys.exit(1)
                else:
                    end_analysis = True
                    break
        if end_analysis:
            break

    logging.debug('Successfully checked for overlapping zones')

    # Group zones by input file
    files_indexes = dict()
    zones_by_file = dict()
    for chrom, chrom_zones in zones.items():
        for zone in chrom_zones:
            for file_index, file in enumerate(zone[2]):
                files_indexes[file] = file_index
                if file not in zones_by_file:
                    zones_by_file[file] = dict()
                if chrom not in zones_by_file[file]:
                    zones_by_file[file][chrom] = []
                zones_by_file[file][chrom].append(zone[:2])

    del zones

    zones_hash = hashlib.md5(str(zones_by_file).encode()).hexdigest()
    output_files_dict = dict()
    pool = ProcessPoolExecutor(args.maximum_processes)
    memory_per_process = args.maximum_memory / args.maximum_processes
    tasks = []
    for filename, zones in zones_by_file.items():
        file_index = files_indexes[filename]
        filename_hash = hashlib.md5(filename.encode() + zones_hash.encode()).hexdigest()
        output_filename = f'{args.outputs[file_index]}_{filename_hash}.{args.outputs[file_index].split(".")[-1]}'
        if file_index not in output_files_dict:
            output_files_dict[file_index] = []
        output_files_dict[file_index].append(output_filename)
        # If output file already exists and does have EOF (not SAM), skip it
        if os.path.exists(output_filename) and not output_filename.endswith('.sam'):
            try:
                pysam.quickcheck(output_filename)
                logging.debug(f'Skipping {filename} because {output_filename} already exists')
                continue
            except:
                pass
        task = pool.submit(_write_zones, zones, file_index, filename,
                           output_filename, args.multithreading, memory_per_process, args.fasta_ref)
        tasks.append(task)
    for task in tasks:
        task.result()
    pool.shutdown()

    # Merge files in parallel
    pool = ProcessPoolExecutor(len(output_files_dict))
    processes_per_file = max(args.maximum_processes // len(output_files_dict), 1)
    for file_index, output_files in output_files_dict.items():
        logging.debug(f'Merging {len(output_files)} files into {args.outputs[file_index]}')
        pool.submit(_merge_files, processes_per_file, output_files, args.outputs[file_index])
    for task in tasks:
        task.result()
    pool.shutdown()

    logging.debug('Done')
