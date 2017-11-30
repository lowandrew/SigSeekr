#!/usr/bin/env python

import multiprocessing
import subprocess
import argparse
import textwrap
import shutil
import glob
import time
import os
import re
from Bio import SeqIO
from biotools import kmc
from biotools import bbtools
from accessoryFunctions.accessoryFunctions import printtime


def find_paired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to try to find paired fastq files. Should be able to find anything fastq.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2.
    :return: List containing pairs of fastq files, in format [[forward_1, reverse_1], [forward_2, reverse_2]], etc.
    """
    pair_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


def find_unpaired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    unpaired_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        if forward_id in name and not os.path.isfile(name.replace(forward_id, reverse_id)):
            unpaired_list.append(name)
        elif forward_id not in name and reverse_id not in name:
            unpaired_list.append(name)
        elif reverse_id in name and not os.path.isfile(name.replace(reverse_id, forward_id)):
            unpaired_list.append(name)
    return unpaired_list


def make_inclusion_kmerdb(inclusion_folder, output_db, forward_id='_R1', reverse_id='_R2', tmpdir='tmpinclusion',
                          maxmem='12'):
    # Make the tmpdir, if it doesn't exist already.
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Get lists of everything - fasta, paired fastq, unpaired fastq.
    fastas = glob.glob(os.path.join(inclusion_folder, '*.f*a'))
    paired_fastqs = find_paired_reads(inclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    unpaired_fastqs = find_unpaired_reads(inclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    # Make a database for each item in each list, and place it into the tmpdir.
    i = 1
    for fasta in fastas:
        kmc.kmc(fasta, os.path.join(tmpdir, 'database{}'.format(str(i))), fm='', m=maxmem)
        i += 1
    for pair in paired_fastqs:
        kmc.kmc(forward_in=pair[0], reverse_in=pair[1], database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2, m=maxmem)  # For fastqs, make min_occurrence two to hopefully filter out sequencing errors.
        i += 1
    for fastq in unpaired_fastqs:
        kmc.kmc(forward_in=fastq, database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2, m=maxmem)  # For fastqs, make min_occurrence two to hopefully filter out sequencing errors.
        i += 1
    # Create a command file to allow kmc to get an intersection of all the inclusion databases created and write to our
    # final inclusion database.
    with open(os.path.join(tmpdir, 'command_file'), 'w') as f:
        f.write('INPUT:\n')
        for j in range(i - 1):
            f.write('set{} = {}\n'.format(str(j + 1), os.path.join(tmpdir, 'database{}'.format(str(j + 1)))))
        f.write('OUTPUT:\n{} = '.format(output_db))
        for j in range(i - 1):
            if j < (i - 2):
                f.write('set{}*'.format(str(j + 1)))
            else:
                f.write('set{}\n'.format(str(j + 1)))
    cmd = 'kmc_tools complex {}'.format(os.path.join(tmpdir, 'command_file'))
    with open('asdf.txt', 'w') as f:
        subprocess.call(cmd, shell=True, stderr=f, stdout=f)
    shutil.rmtree(tmpdir)


def make_exclusion_kmerdb(exclusion_folder, output_db, forward_id='_R1', reverse_id='_R2', tmpdir='tmpexclusion',
                          maxmem='12'):
    # Make the tmpdir, if it doesn't exist already.
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Get lists of everything - fasta, paired fastq, unpaired fastq.
    fastas = glob.glob(os.path.join(exclusion_folder, '*.f*a'))
    paired_fastqs = find_paired_reads(exclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    unpaired_fastqs = find_unpaired_reads(exclusion_folder, forward_id=forward_id, reverse_id=reverse_id)
    # Make a database for each item in each list, and place it into the tmpdir.
    i = 1
    for fasta in fastas:
        kmc.kmc(fasta, os.path.join(tmpdir, 'database{}'.format(str(i))), fm='', m=maxmem)
        i += 1
    for pair in paired_fastqs:
        kmc.kmc(forward_in=pair[0], reverse_in=pair[1], database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2, m=maxmem)  # For fastqs, make min_occurrence two to hopefully filter out sequencing errors.
        i += 1
    for fastq in unpaired_fastqs:
        kmc.kmc(forward_in=fastq, database_name=os.path.join(tmpdir, 'database{}'.format(str(i))),
                min_occurrences=2, m=maxmem)  # For fastqs, make min_occurrence two to hopefully filter out sequencing errors.
        i += 1
    # Create a command file to allow kmc to do a union of all the databases you've created and write them to our final
    # exclusion db.
    with open(os.path.join(tmpdir, 'command_file'), 'w') as f:
        f.write('INPUT:\n')
        for j in range(i - 1):
            f.write('set{} = {}\n'.format(str(j + 1), os.path.join(tmpdir, 'database{}'.format(str(j + 1)))))
        f.write('OUTPUT:\n{} = '.format(output_db))
        for j in range(i - 1):
            if j < (i - 2):
                f.write('set{}+'.format(str(j + 1)))
            else:
                f.write('set{}\n'.format(str(j + 1)))
    cmd = 'kmc_tools complex {}'.format(os.path.join(tmpdir, 'command_file'))
    with open('asdf.txt', 'w') as f:
        subprocess.call(cmd, shell=True, stderr=f, stdout=f)
    shutil.rmtree(tmpdir)


def kmers_to_fasta(kmer_file, output_fasta):
    with open(kmer_file) as infile:
        lines = infile.readlines()
    with open(output_fasta, 'w') as outfile:
        i = 1
        for line in lines:
            sequence = line.split()[0]  # Sequence is the first thing in the split
            outfile.write('>kmer{}\n'.format(str(i)))
            outfile.write(sequence + '\n')
            i += 1


def remove_n(input_fasta, output_fasta):
    contigs = SeqIO.parse(input_fasta, 'fasta')
    j = 1
    for contig in contigs:
        sequence = str(contig.seq)
        uniques = re.split('N+', sequence)
        with open(output_fasta, 'a+') as outfile:
            i = 1
            for unique in uniques:
                if unique != '':
                    outfile.write('>contig{}_sequence{}\n'.format(str(j), str(i)))
                    unique = textwrap.fill(unique)
                    outfile.write(unique + '\n')
                    i += 1
            j += 1


def replace_by_index(stretch, seq):
    stretch = stretch.split(':')
    start = int(stretch[0])
    end = int(stretch[1])
    seq = seq[:start] + 'N'*(end-start) + seq[end:]
    return seq


def mask_fasta(input_fasta, output_fasta, bedfile):
    to_mask = dict()
    with open(bedfile) as bed:
        lines = bed.readlines()
    for line in lines:
        line = line.rstrip()
        x = line.split()
        coverage = x[-1]
        end = x[-2]
        start = x[-3]
        name = ' '.join(x[:-3])
        if coverage == '0':
            if name in to_mask:
                to_mask[name].append(start + ':' + end)
            else:
                to_mask[name] = [start + ':' + end]
    fasta_in = SeqIO.parse(input_fasta, 'fasta')
    for contig in fasta_in:
        seq = str(contig.seq)
        if contig.description in to_mask:
            for item in to_mask[contig.description]:
                seq = replace_by_index(item, seq)
            with open(output_fasta, 'a+') as outfile:
                outfile.write('>{}\n'.format(contig.description))
                outfile.write(seq + '\n')


def generate_bedfile(ref_fasta, kmers, output_bedfile, tmpdir='bedgentmp'):
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # First, need to generate a bam file - align the kmers to a reference fasta genome.
    bbtools.bbmap(ref_fasta, kmers, os.path.join(tmpdir, 'out.bam'))
    # Once the bam file is generated, turn it into a sorted bamfile so that bedtools can work with it.
    cmd = 'samtools sort {bamfile} -o {sorted_bamfile}'.format(bamfile=os.path.join(tmpdir, 'out.bam'),
                                                               sorted_bamfile=os.path.join(tmpdir, 'out_sorted.bam'))
    subprocess.call(cmd, shell=True)
    # Use bedtools to get genome coverage, so that we know what to mask.
    cmd = 'bedtools genomecov -ibam {sorted_bamfile} -bga' \
          ' > {output_bed}'.format(sorted_bamfile=os.path.join(tmpdir, 'out_sorted.bam'),
                                   output_bed=output_bedfile)
    subprocess.call(cmd, shell=True)
    shutil.rmtree(tmpdir)


def main(args):
    start = time.time()
    # Make the necessary inclusion and exclusion kmer sets.
    printtime('Creating inclusion kmer set...', start)
    make_inclusion_kmerdb(args.inclusion, os.path.join(args.output_folder, 'inclusion_db'))
    printtime('Creating exclusion kmer set...', start)
    make_exclusion_kmerdb(args.exclusion, os.path.join(args.output_folder, 'exclusion_db'))
    # Now start trying to subtract kmer sets, see how it goes.
    exclusion_cutoff = 1
    while exclusion_cutoff < 10:
        printtime('Subtracting exclusion kmers from inclusion kmers with cutoff {}...'.format(str(exclusion_cutoff)), start)
        kmc.subtract(os.path.join(args.output_folder, 'inclusion_db'), os.path.join(args.output_folder, 'exclusion_db'),
                     os.path.join(args.output_folder, 'unique_to_inclusion_db'), exclude_below=exclusion_cutoff)
        kmc.dump(os.path.join(args.output_folder, 'unique_to_inclusion_db'),
                 os.path.join(args.output_folder, 'unique_kmers.txt'))
        # Now need to check if any kmers are present, and if not, increment the counter to allow a more lax search.
        with open(os.path.join(args.output_folder, 'unique_kmers.txt')) as f:
            lines = f.readlines()
        if lines != []:
            break
        exclusion_cutoff += 1
    # Convert our kmers to FASTA format for usage with other programs.
    kmers_to_fasta(os.path.join(args.output_folder, 'unique_kmers.txt'),
                   os.path.join(args.output_folder, 'inclusion_kmers.fasta'))
    # Now that we have kmers that are unique to inclusion sequence, need to map them back to an inclusion genome, if
    # we have one in FASTA format. This will allow us to find unique regions, instead of just kmers.
    if len(glob.glob(os.path.join(args.inclusion, '*.f*a'))) > 0:
        ref_fasta = glob.glob(os.path.join(args.inclusion, '*.f*a'))[0]
        # Get inclusion kmers into FASTA format.
        generate_bedfile(ref_fasta, os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                         os.path.join(args.output_folder, 'regions_to_mask.bed'))
        mask_fasta(ref_fasta, os.path.join(args.output_folder, 'inclusion_sequence.fasta'),
                   os.path.join(args.output_folder, 'regions_to_mask.bed'))
        remove_n(os.path.join(args.output_folder, 'inclusion_sequence.fasta'),
                 os.path.join(args.output_folder, 'sigseekr_result.fasta'))
    # If we want to find PCR primers, try to filter out any inclusion kmers that are close to exclusion kmers.
    if args.pcr:
        # First step is to create fasta of exclusion kmers.
        kmc.dump(os.path.join(args.output_folder, 'exclusion_db'),
                 os.path.join(args.output_folder, 'exclusion_kmers.txt'))
        kmers_to_fasta(os.path.join(args.output_folder, 'exclusion_kmers.txt'),
                       os.path.join(args.output_folder, 'exclusion_kmers.fasta'))
        # Now use bbduk with small kmer size (k=15) to filter out inclusion kmers that have exclusions that are close.
        bbtools.bbduk_filter(reference=os.path.join(args.output_folder, 'exclusion_kmers.fasta'),
                             forward_in=os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                             forward_out=os.path.join(args.output_folder, 'pcr_kmers.fasta'),
                             k='15')
        # Next step: Get distances between potential primers by mapping back to a reference (if it exists) and getting
        # distances.


if __name__ == '__main__':
    start = time.time()
    num_cpus = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inclusion',
                        type=str,
                        required=True,
                        help='Path to folder containing genome(s) you want signature sequences for.'
                             ' Genomes must be in FASTA format (.fasta), and should not be compressed.')
    parser.add_argument('-e', '--exclusion',
                        type=str,
                        required=True,
                        help='Path to folder containing exclusion genome(s) - those you do not want signature'
                             ' sequences for. Genomes must be in FASTA format (.fasta) and should not be compressed.')
    parser.add_argument('-o', '--output_folder',
                        type=str,
                        required=True,
                        help='Path to folder where you want to store output files. Folder will be created if it '
                             'does not exist.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=num_cpus,
                        help='Number of threads to run analysis on. Defaults to number of cores on your machine.')
    parser.add_argument('-pcr', '--pcr',
                        default=False,
                        action='store_true',
                        help='Enable to filter out inclusion kmers that have close relatives in exclusion kmers.')
    args = parser.parse_args()
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    main(args)
