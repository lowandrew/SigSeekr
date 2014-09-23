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
from collections import defaultdict
import os, glob, GeneSeekr, shutil, json, time

def retriever(genomes, output):
    if not os.path.exists(output + "Genomes"):
        os.mkdir(output + "Genomes")
    for folders in glob.glob(genomes + "/*"):
        if os.path.exists(folders + "/Best_Assemblies"):
            for fasta in glob.glob(folders + "/Best_Assemblies/*"):
                shutil.copy(fasta, output + "Genomes")

def jsonUpGoer(jsonfile, markers, genomes, outdir):
    if os.path.isfile(jsonfile):
        genedict = json.load(open(jsonfile))
    else:
        genedict = GeneSeekr.blaster(markers, genomes, outdir, "USSpip")
        json.dump(genedict, open(jsonfile, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    return genedict

def compareType(TargetrMLST, nonTargetrMLST):
    '''
    Compare the rMLST types by making sure that the alleles are exactly the same size
    Record list of genomes not included in non-target selection
    Requires target and non-target dictionaries as inputs
    '''
    typing = defaultdict(dict)
    removed = defaultdict(list)
    for genome in TargetrMLST:
        if genome not in typing:
            typing[genome] = defaultdict(int)
        for gene in sorted(TargetrMLST[genome]):
            for nontarget in nonTargetrMLST:
                if nontarget not in typing[genome]:  # if nontarget genome not in typing dictionary then add it
                    typing[genome][nontarget] = 0
                if TargetrMLST[genome][gene] == nonTargetrMLST[genome][gene]:
                    typing[genome][nontarget] += 1
    for genome in typing:
        for nontarget in typing[genome]:
            if typing[genome][nontarget] == 53: # Actual number of alleles ... I hate you Keith
                removed[genome].append(nontarget)
                typing[genome].pop(nontarget)
    return typing, removed




def sorter(markers, genomes, outdir, target):
    '''Strip first allele off each locus to feed into geneseekr and return dictionary
    '''
    smallMLST = outdir + "rMLST/"
    if not os.path.exists(outdir + "Genomes/"):
        retriever(genomes, outdir)
    if not os.path.exists(outdir + "tmp/"):
        os.mkdir(outdir + "tmp/")
    genomes = outdir + "Genomes/"
    jsonfile = '%sgenedict.json' %  outdir
    nonTargetrMLST = jsonUpGoer(jsonfile, markers, genomes, outdir)
    if os.path.isdir(target):  # Determine if target is a folder
        targets = glob.glob(target + "*")
        targetjson = '%stargetdict.json' % target
    elif os.path.isfile(target):
        targets = target
        targetjson = '%stargetdict.json' % outdir
    else:
        print "The variable \"--targets\" is not a folder or file "
        return
    # TargetrMLST = jsonUpGoer(targetjson, targets, genomes, outdir)
    # typing, removed = compareType(TargetrMLST, nonTargetrMLST)
    # json.dump(typing, open(outdir + 'typing.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    # json.dump(removed, open(outdir + 'removed.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))








#Parser for arguments
parser = ArgumentParser(description='Find Universal Strain-Specifc Probes')
parser.add_argument('--version', action='version', version='%(prog)s v0.5')
parser.add_argument('-o', '--output', required=True, help='Specify output directory')
parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify rMLST folder')
parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
# parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
args = vars(parser.parse_args())

sorter(args['marker'], args['input'], args['output'], args['target'])