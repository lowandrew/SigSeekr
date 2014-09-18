__author__ = 'mikeknowles'
""" This is the core script, its mission:
Retrieve genomes
This will require the user to download rMLST data
To sort rMLST results and remove closely related sequences
then to prepare the data for strain-specific probe idenification
"""
from argparse import ArgumentParser
from Bio import SeqIO
from textwrap import fill
import os, glob, GeneSeekr, shutil, json

def retriever(genomes, output):
    if not os.path.exists(output + "Genomes"):
        os.mkdir(output + "Genomes")
    for folders in glob.glob(genomes + "/*"):
        if os.path.exists(folders + "/Best_Assemblies"):
            for fasta in glob.glob(folders + "/Best_Assemblies/*"):
                shutil.copy(fasta, output + "Genomes")


def sorter(markers, genomes, outdir):
    '''Strip first allele off each locus to feed into geneseekr and return dictionary
    '''
    smallMLST = outdir + "rMLST/"
    if not os.path.exists(outdir + "Genomes/"):
        retriever(genomes, outdir)
    if not os.path.exists(outdir + "tmp/"):
        os.mkdir(outdir + "tmp/")
    genomes = outdir + "Genomes/"

    # if not os.path.exists(smallMLST):
    #     os.mkdir(smallMLST)
    # for fasta in glob.glob(markers + "/*.fas"):
    #     for record in SeqIO.parse(fasta, "fasta"):
    #         with open("%s%s.fasta" %(smallMLST, record.id[0:10]), 'w') as handle:
    #             handle.write(">%s\n%s" % (record.id, fill(str(record.seq), width=80)))
    #         break
    jsonfile = '%sgenedict.json' %  markers
    if os.path.isfile(jsonfile):
        genedict = json.load(open(jsonfile))
    else:
        genedict = GeneSeekr.blaster(markers, genomes, outdir, 0, "USSpip")
        json.dump(genedict, open(jsonfile, 'w'), sort_keys=True, indent=4, separators=(',', ': '))




#Parser for arguments
parser = ArgumentParser(description='Find Universal Strain-Specifc Probes')
parser.add_argument('--version', action='version', version='%(prog)s v0.5')
parser.add_argument('-o', '--output', required=True, help='Specify output directory')
parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify rMLST folder')
args = vars(parser.parse_args())

sorter(args['marker'], args['input'], args['output'])