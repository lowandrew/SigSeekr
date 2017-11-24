#!/usr/bin/env python

import subprocess
import argparse
import textwrap
import glob
import time
import os
import re
from Bio import SeqIO
from biotools import kmc
from biotools import bbtools
from accessoryFunctions.accessoryFunctions import printtime


"""
Step 1: Concatenate inclusion and exclusion genomes each into a big multifasta.
Step 2: Kmerize each multifasta.
Step 3: Subtract kmers in exclusion kmer database from kmers in inclusion database.
Step 4: Do assembly on kmers unique to inclusion strains to get unique regions.

May have to iterate on kmer size: Start long, and then move towards smaller kmers if it ends up failing.
"""


def find_inclusion_kmers(inclusion_folder, inclusion_db):
    inclusion_fastas = glob.glob(os.path.join(inclusion_folder, '*fasta'))
    if len(inclusion_fastas) == 1:
        kmc.kmc(inclusion_fastas[0], inclusion_db, fm='')
    elif len(inclusion_fastas) == 2:
        # Will need to get rid of these files at some point.
        kmc.kmc(inclusion_fastas[0], os.path.join(inclusion_folder, 'database1'), fm='')
        kmc.kmc(inclusion_fastas[1], os.path.join(inclusion_folder, 'database2'), fm='')
        kmc.intersect(os.path.join(inclusion_folder, 'database1'),
                      os.path.join(inclusion_folder, 'database2'),
                      inclusion_db)
    elif len(inclusion_fastas) > 2:
        # Create a database for each fasta.
        i = 1
        for fasta in inclusion_fastas:
            kmc.kmc(fasta, os.path.join(inclusion_folder, 'database' + str(i)), fm='')
            i += 1
        # Now need to do an intersection of every database created.
        # Step 1: Creaate a KMC command file.
        with open(os.path.join(inclusion_folder, 'command_file'), 'w') as f:
            f.write('INPUT:\n')
            for i in range(len(inclusion_fastas)):
                f.write('set{} = {}\n'.format(str(i + 1), os.path.join(inclusion_folder, 'database' + str(i + 1))))
            f.write('OUTPUT:\n{} = '.format(inclusion_db))
            for i in range(len(inclusion_fastas)):
                if i < len(inclusion_fastas) - 1:
                    f.write('set{}*'.format(str(i + 1)))
                else:
                    f.write('set{}\n'.format(str(i + 1)))
        cmd = 'kmc_tools complex {}'.format(os.path.join(inclusion_folder, 'command_file'))
        subprocess.call(cmd, shell=True)
        os.remove(os.path.join(inclusion_folder, 'command_file'))
    # Do cleanup on any file that may have been created.
    cleanup_files = glob.glob(os.path.join(inclusion_folder, 'database*kmc*'))
    for name in cleanup_files:
        os.remove(name)


def concatenate_fastas(fasta_folder, output_fasta):
    # This is slow - see if you can make if faster.
    fasta_files = glob.glob(os.path.join(fasta_folder, '*fasta'))
    with open(output_fasta, 'w') as outfile:
        for fasta in fasta_files:
            with open(fasta) as infile:
                outfile.write(infile.read())


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
    # Need to know how long strings of Ns between unique regions are - will attempt to find out exactly how to do this
    # tomorrow.
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


if __name__ == '__main__':
    start = time.time()
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
    args = parser.parse_args()
    # Make output folder if it doesn't already exist.
    if not os.path.isdir(args.output_folder):
        os.makedirs(args.output_folder)
    # Concatenate exclusion sequences.
    printtime('Finding kmers common to inclusion sequences', start)
    find_inclusion_kmers(args.inclusion, os.path.join(args.output_folder, 'inclusion_db'))
    # Concatenate exclusion kmers.
    printtime('Concatenating exclusion sequences', start)
    concatenate_fastas(args.exclusion, os.path.join(args.output_folder, 'exclusion.fasta'))
    # Run KMC on inclusion and exclusion fastas.
    # kmc.kmc(os.path.join(args.output_folder, 'inclusion.fasta'),
    #         os.path.join(args.output_folder, 'inclusion_db'), fm='')
    printtime('Finding exclusion kmers', start)
    kmc.kmc(os.path.join(args.output_folder, 'exclusion.fasta'),
            os.path.join(args.output_folder, 'exclusion_db'), fm='')
    # Subtract exclusion kmers from inclusion kmers to find only kmers that are unique to inclusion.
    printtime('Determining which kmers are unique to inclusion sequences', start)
    kmc.subtract(os.path.join(args.output_folder, 'inclusion_db'),
                 os.path.join(args.output_folder, 'exclusion_db'),
                 os.path.join(args.output_folder, 'unique_to_inclusion'))
    # Dump unique to inclusion kmers into a file.
    printtime('Extracting Signature Sequences...', start)
    kmc.dump(os.path.join(args.output_folder, 'unique_to_inclusion'),
             os.path.join(args.output_folder, 'inclusion_kmers.txt'))
    # Convert kmers file to a fasta file.
    kmers_to_fasta(os.path.join(args.output_folder, 'inclusion_kmers.txt'),
                   os.path.join(args.output_folder, 'inclusion_kmers.fasta'))
    # Map kmers back to the inclusion fasta
    fasta_to_use = glob.glob(os.path.join(args.inclusion, '*fasta'))[0]
    bbtools.bbmap(os.path.join(fasta_to_use),
                  os.path.join(args.output_folder, 'inclusion_kmers.fasta'),
                  os.path.join(args.output_folder, 'out.bam'))
    # Sort the resulting bam file.
    cmd = 'samtools sort {bamfile} -o {sorted_bamfile}'.format(bamfile=os.path.join(args.output_folder, 'out.bam'),
                                                               sorted_bamfile=os.path.join(args.output_folder, 'out_sorted.bam'))
    subprocess.call(cmd, shell=True)
    # Use bedtools to get genome coverage for each region of inclusion file.
    cmd = 'bedtools genomecov -ibam {sorted_bamfile} -bga' \
          ' > {output_bed}'.format(sorted_bamfile=os.path.join(args.output_folder, 'out_sorted.bam'),
                                   output_bed=os.path.join(args.output_folder, 'bedfile.bed'))
    subprocess.call(cmd, shell=True)
    # Mask the regions of the fasta file that have zero coverage.
    mask_fasta(fasta_to_use, os.path.join(args.output_folder, 'masked_sequence.fasta'),
               os.path.join(args.output_folder, 'bedfile.bed'))
    # Replace all Ns with nothing, so that we're left with only unique sequence.
    remove_n(os.path.join(args.output_folder, 'masked_sequence.fasta'),
             os.path.join(args.output_folder, 'unique_sequence.fasta'))
    printtime('Signature sequence finding complete :D', start, '\033[0;92m')
