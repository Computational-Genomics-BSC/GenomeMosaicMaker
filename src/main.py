# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License

import os
import sys
import argparse
import bisect
import hashlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import logging
import pysam


from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantType

DEFAULT_RG_NAME = 'COMBINED'


def _write_read(output_file_obj, read, input_filename, split_rg):
    new_query_name = hashlib.md5((f'{read.query_name}_{input_filename}').encode()).hexdigest()
    read.query_name = new_query_name
    if split_rg:
        # Set read group if not present
        rg_tag = read.get_tag('RG')
        if rg_tag is None:
            rg_tag = DEFAULT_RG_NAME
        read.set_tag('RG', hashlib.md5((rg_tag).encode()).hexdigest())
    else:
        read.set_tag('RG', DEFAULT_RG_NAME)
    output_file_obj.write(read)


def _get_reads_in_zones_chunk(input_file, zone_chunk, reference_fasta):
    input_file_obj = pysam.AlignmentFile(input_file, reference_filename=reference_fasta)
    inzone_reads = set()
    for chrom, start, end in zone_chunk:
        for read in input_file_obj.fetch(chrom, start, end):
            inzone_reads.add(read.query_name)
    return inzone_reads


def _get_reads_in_zones(input_file, zones_dict, num_processes, reference_fasta):
    # Get all zones in a list
    zone_list = []
    for chrom, chrom_zones in zones_dict.items():
        for zone in chrom_zones:
            zone_list.append((chrom, zone[0], zone[1]))
    # Separate zones in chunks
    zone_chunks = []
    chunk_size = int(len(zone_list) / num_processes) + 1
    for i in range(0, len(zone_list), chunk_size):
        zone_chunks.append(zone_list[i:i + chunk_size])
    # Get reads in zones
    inzone_reads = set()
    with ThreadPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for zone_chunk in zone_chunks:
            futures.append(executor.submit(_get_reads_in_zones_chunk, input_file, zone_chunk, reference_fasta))
        for future in futures:
            inzone_reads.update(future.result())
    return inzone_reads


def _get_reads_generator(input_file_objs, ref_names_order):
    # Generator that yields reads from all input files in order
    # Get all files read generators
    read_generators = []
    for input_file_obj in input_file_objs:
        read_generators.append(input_file_obj.fetch(until_eof=True))
    reads = [next(read_generator) for read_generator in read_generators]
    # Yield reads in order (by reference name and then by start position)
    while reads:
        if len(reads) == 1:
            yield reads[0]
            try:
                reads[0] = next(read_generators[0])
            except StopIteration:
                reads.pop(0)
                read_generators.pop(0)
            continue
        # Get lowest reference name in the reference names order in the read list
        min_ref_name = min(reads, key=lambda read: ref_names_order[read.reference_name]).reference_name
        # Get the index of the read with the lowest start position in with the lowest reference name
        min_start_index = min(
            range(len(reads)), key=lambda i: reads[i].reference_start if reads[i].reference_name == min_ref_name else sys.maxsize)
        # Yield the read
        yield reads[min_start_index]
        # Get the next read
        try:
            reads[min_start_index] = next(read_generators[min_start_index])
        except StopIteration:
            reads.pop(min_start_index)
            read_generators.pop(min_start_index)
    yield None


def _write_zones(zones, input_filename, output_filename, available_threads, fasta_ref=None):
    # If .done file exists, skip
    if os.path.isfile(f'{output_filename}.done'):
        logging.info(f'Skipping {input_filename} > {output_filename}')
        return
    # Get reads in zones
    qnames_set = _get_reads_in_zones(input_filename, zones, available_threads, fasta_ref)
    del zones

    # Read the whole file and write the reads in the zones
    with pysam.AlignmentFile(input_filename, 'r', threads=available_threads, reference_filename=fasta_ref) as input_file_obj:
        with _open_output_file(output_filename, [input_file_obj.header.to_dict()], available_threads, fasta_ref) as output_file_obj:
            for read in input_file_obj.fetch(until_eof=True):
                if read.query_name in qnames_set:
                    output_file_obj.write(read)
    logging.debug(f'Loaded reads for {input_filename} > {output_filename}')
    # Write .done file to indicate that the file has been processed
    with open(f'{output_filename}.done', 'w') as done_file:
        done_file.write('')


def _get_sam_header(file_obj, file_index, split_rg, rewrite_rg=True):
    header = file_obj.header.to_dict()
    # Remove the PG lines
    header['PG'] = []
    old_read_group_dict = {}
    for rg in header.get('RG', []):
        old_read_group_dict[rg['ID']] = rg
    if rewrite_rg:
        if split_rg:
            # Set the RG line
            new_read_groups = []
            # Add the default RG line
            old_read_group_dict[DEFAULT_RG_NAME] = {'ID': DEFAULT_RG_NAME}
            # Add the SM tag to each RG and remove the other entries
            for read_group in old_read_group_dict.values():
                new_read_groups.append(
                    {'ID': hashlib.md5((read_group['ID']).encode()).hexdigest(), 'SM': str(file_index)})
        else:
            new_read_groups = [{'ID': DEFAULT_RG_NAME, 'SM': str(file_index)}]
    else:
        new_read_groups = list(old_read_group_dict.values())
    header['RG'] = new_read_groups
    return header


def _open_output_file(output_file, headers, num_threads, fasta_ref):
    # Get write mode based on the output file extension
    if output_file.endswith('.bam'):
        write_mode = 'wb'
    elif output_file.endswith('.sam'):
        write_mode = 'wh'
    elif output_file.endswith('.cram'):
        write_mode = 'wc'
    else:
        raise Exception('Invalid output file extension')

    # Combine the RG from the headers
    new_header = headers[0]
    read_group_dict = {}
    for header in headers:
        for rg in header.get('RG', []):
            read_group_dict[rg['ID']] = rg
    new_header['RG'] = list(read_group_dict.values())
    # Open the output file
    return pysam.AlignmentFile(output_file, write_mode, header=new_header, threads=num_threads, reference_filename=fasta_ref)


def _merge_files(temp_output_files, file_index, output_file, zones, num_processes, fasta_ref, split_rg, canvas_file=None):
    temp_output_files_objs = [pysam.AlignmentFile(f, reference_filename=fasta_ref, threads=num_processes)
                              for f in temp_output_files]
    template_headers = [_get_sam_header(f_obj, file_index, split_rg, rewrite_rg=True) for f_obj in temp_output_files_objs] + \
        ([_get_sam_header(pysam.AlignmentFile(canvas_file, reference_filename=fasta_ref), file_index, split_rg, rewrite_rg=True)]
         if canvas_file is not None else [])
    # Open the output file
    output_file_obj = _open_output_file(output_file, template_headers, num_processes, fasta_ref)
    # Get reference names order
    ref_names_order = dict()
    for idx, ref_name in enumerate(output_file_obj.references):
        ref_names_order[ref_name] = idx
    ref_names_order['*'] = len(ref_names_order)
    # Get generator of reads from the temp output files
    output_reads_generator = _get_reads_generator(temp_output_files_objs, ref_names_order)
    # Get the first read from the generator
    output_read = next(output_reads_generator)
    if canvas_file is not None:
        # Get the reads from the canvas file that are in the zones
        inzone_reads = _get_reads_in_zones(canvas_file, zones, num_processes, fasta_ref)
        # Open the canvas file
        canvas_file_obj = pysam.AlignmentFile(canvas_file, threads=num_processes, reference_filename=fasta_ref)
        # Write the reads from the canvas file
        for canvas_read in canvas_file_obj.fetch(until_eof=True):
            # Skip reads that are not in the original reference
            if canvas_read.reference_name not in ref_names_order or canvas_read.next_reference_name not in ref_names_order:
                continue
            # Write the output reads until we reach the canvas read
            while output_read is not None and \
                (ref_names_order[output_read.reference_name] < ref_names_order[canvas_read.reference_name] or
                 (output_read.reference_name == canvas_read.reference_name and output_read.reference_start < canvas_read.reference_start)):
                _write_read(output_file_obj, output_read, canvas_file, split_rg)
                output_read = next(output_reads_generator)
            # Avoid writing the canvas read if it is in a zone
            if canvas_read.query_name in inzone_reads:
                continue
            # Write the canvas read
            _write_read(output_file_obj, canvas_read, canvas_file, split_rg)
        # Close the canvas file
        canvas_file_obj.close()
    # Write the remaining output reads
    while output_read is not None:
        _write_read(output_file_obj, output_read, canvas_file, split_rg)
        output_read = next(output_reads_generator)

    # Close the output file
    output_file_obj.close()
    # Close the temp output files
    _ = [f.close() for f in temp_output_files_objs]


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
        if var_type == VariantType.SNV.name:
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files, str(var_obj)))
        elif var_type == VariantType.INS.name:
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] +
                          row['length'] + padding, files, str(var_obj)))
        elif var_type == VariantType.TRA.name:
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files, str(var_obj)))
            bisect.insort(zones[end_chrom], (row['end'] - padding, row['end'] + padding, files, str(var_obj)))
        elif var_type == VariantType.INV.name:
            if row['start'] + padding > row['end'] - padding:
                bisect.insort(zones[start_chrom], (row['start'] - padding, row['end'] + padding, files, str(var_obj)))
            else:
                bisect.insort(zones[start_chrom], (row['start'] - padding, row['start'] + padding, files, str(var_obj)))
                bisect.insort(zones[end_chrom], (row['end'] - padding, row['end'] + padding, files, str(var_obj)))
        else:
            bisect.insort(zones[start_chrom], (row['start'] - padding, row['end'] + padding, files, str(var_obj)))
    return zones


def write_zones_bed_file(zones, output_bed, extra_padding=0):
    output_bed = args.outputs[0] + '_zones.bed'
    with open(output_bed, 'w') as bed_file:
        for chrom, chrom_zones in zones.items():
            for zone in chrom_zones:
                bed_file.write(f'{chrom}\t{max(0, zone[0] - extra_padding)}\t{zone[1] + extra_padding}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, help='Input VCF file')
    parser.add_argument('--outputs', '-o', required=True, nargs='+', help='Output alignment files')
    parser.add_argument('--canvas-files', '-if', nargs='+', help='canvas alignment files')
    parser.add_argument('--padding', '-p', type=int, default=1000, help='Padding around the variants')
    parser.add_argument('--maximum-processes', '-mp', type=int, default=1,
                        help='Maximum number of physical processes to use')
    parser.add_argument('--fasta-ref', '-f', type=str, help='Fasta reference file (used for CRAM files)')
    parser.add_argument('--split-read-groups', action='store_true',
                        help='Keep the read groups separate (they will be anonymized)')

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG)

    # Convert everything to absolute paths
    args.input = os.path.abspath(args.input)
    args.outputs = [os.path.abspath(output) for output in args.outputs]

    if args.canvas_files is not None:
        args.canvas_files = [os.path.abspath(canvas_file) for canvas_file in args.canvas_files]
        # Make sure the same number of canvas files and output files are provided
        if len(args.canvas_files) != len(args.outputs):
            raise Exception('The same number of canvas files and output files must be provided')
        # Make sure the canvas files and output files are not the same
        for canvas_file, output_file in zip(args.canvas_files, args.outputs):
            if canvas_file == output_file:
                raise Exception('The canvas files and output files must be different')
        # Make sure the canvas files exist and are not empty
        for canvas_file in args.canvas_files:
            if not os.path.isfile(canvas_file) or os.path.getsize(canvas_file) == 0:
                raise Exception(f'The canvas file {canvas_file} does not exist or is empty')

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

    # Create a BED file with the zones
    output_bed = args.outputs[0] + '_zones.bed'
    write_zones_bed_file(zones, output_bed, args.padding)
    logging.debug(f'Wrote zones to {output_bed}')

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

    # Write zones to files
    output_files_dict = dict()
    pool = ProcessPoolExecutor(min(args.maximum_processes, len(zones_by_file)))
    processes_per_file = max(args.maximum_processes // len(zones_by_file), 1)
    tasks = []
    for filename, file_zones in zones_by_file.items():
        zones_hash = hashlib.md5(str(file_zones).encode()).hexdigest()
        file_index = files_indexes[filename]
        filename_hash = hashlib.md5(filename.encode() + zones_hash.encode()).hexdigest()
        output_filename = f'{args.outputs[file_index]}_{filename_hash}.{args.outputs[file_index].split(".")[-1]}'
        if file_index not in output_files_dict:
            output_files_dict[file_index] = []
        output_files_dict[file_index].append(output_filename)
        task = pool.submit(_write_zones, file_zones, filename, output_filename, processes_per_file, args.fasta_ref)
        tasks.append(task)
    for task in tasks:
        task.result()

    # Merge files in parallel
    processes_per_file = max(args.maximum_processes // len(output_files_dict), 1)
    tasks = []
    for file_index, temp_output_files in output_files_dict.items():
        logging.debug(f'Merging {len(temp_output_files)} files into {args.outputs[file_index]}')
        canvas_file = None if args.canvas_files is None else args.canvas_files[file_index]
        task = pool.submit(_merge_files, temp_output_files, file_index,
                           args.outputs[file_index], zones, processes_per_file, args.fasta_ref, args.split_read_groups, canvas_file)
        tasks.append(task)
    for task in tasks:
        task.result()
    pool.shutdown()

    # Remove temporary files
    for temp_output_files in output_files_dict.values():
        for temp_output_file in temp_output_files:
            os.remove(temp_output_file)
            # Remove .done file
            os.remove(f'{temp_output_file}.done')

    logging.debug('Done')
