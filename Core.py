__author__ = 'mikeknowles'
""" This is the core script, its mission:
Retrieve genomes
This will require the user to download rMLST data
To sort rMLST results and remove closely related sequences
then to prepare the data for strain-specific probe idenification
"""
from argparse import ArgumentParser
from copy import deepcopy
from collections import defaultdict
from glob import glob
from USSPpip import SigSeekr
import os, GeneSeekr, shutil, json

def retriever(genomes, output):
    if not os.path.exists(output + "Genomes"):
        os.mkdir(output + "Genomes")
    for folders in glob(genomes + "/*"): 
        if os.path.exists(folders + "/Best_Assemblies"):
            for fasta in glob(folders + "/Best_Assemblies/*"):
                shutil.copy(fasta, output + "Genomes")

def jsonUpGoer(jsonfile, markers, genomes, outdir, method):
    if os.path.isfile(jsonfile):
        genedict = json.load(open(jsonfile))
    else:
        genedict = GeneSeekr.blaster(markers, genomes, outdir, "USSpip_" + method)
        handle = open(jsonfile, 'w')
        json.dump(genedict, handle, sort_keys=True, indent=4, separators=(',', ': '))
        handle.close()
    return genedict

def alleledictbuild(TargetrMLST, nonTargetrMLST):
    typing = defaultdict(dict)
    for genome in TargetrMLST:  # genome refers to target genome
        if genome not in typing:  # add the target genome to the dictionary
            typing[genome] = defaultdict(int)

        for gene in sorted(TargetrMLST[genome]):  # gene is the rST
            for nontarget in nonTargetrMLST:  # check against this genome
                if nontarget not in typing[genome]:  # if nontarget genome not in typing dictionary then add it
                    typing[genome][nontarget] = 0
                if gene in TargetrMLST[genome] and gene in nonTargetrMLST[nontarget]:  #
                    match = 0
                    for allele in TargetrMLST[genome][gene]:  # multiple allele types possibly present
                        if allele in nonTargetrMLST[nontarget][gene]:
                            match += 1
                    if match == len(nonTargetrMLST[nontarget][gene]):
                        typing[genome][nontarget] += 1
    return typing

def compareType(TargetrMLST, nonTargetrMLST):
    '''
    Compare the rMLST types by making sure that the alleles are exactly the same size
    Record list of genomes not included in non-target selection
    Requires target and non-target dictionaries as inputs
    '''
    removed = defaultdict(list)
    # nontyping = alleledictbuild(nonTargetrMLST, nonTargetrMLST)
    typing = alleledictbuild(TargetrMLST, nonTargetrMLST)
    typing_bak = deepcopy(typing)
    typing = defaultdict(list)
    for genome in typing_bak:
        for nontarget, value in sorted(typing_bak[genome].iteritems(), key=lambda (k, v): (v, k), reverse=True):
        # sort dic*.fionary by alleles in common
            if 0 < value < 53:
                if genome not in typing:
                    typing[genome] = []
                # Actual number of alleles: 53
                typing[genome].append((nontarget, value))
            else:
                removed[genome].append((nontarget, value))
    return typing, removed

def sorter(markers, genomes, outdir, target, evalue, estop):
    '''Strip first allele off each locus to feed into geneseekr and return dictionary
    '''
    # if not os.path.exists(outdir + "Genomes/"):
    #     retriever(genomes, outdir)
    if not os.path.exists(outdir + "tmp/"):
        os.mkdir(outdir + "tmp/")
    # genomes = outdir + "Genomes/"
    nontargets = glob(genomes + "*.fa*")
    jsonfile = '%sgenedict.json' % outdir
    nonTargetrMLST = jsonUpGoer(jsonfile, markers, genomes, outdir, 'nontarget')
    if os.path.exists(target):  # Determine if target is a folder
        targetjson = '%stargetdict.json' % outdir
    else:
        print "The variable \"--targets\" is not a folder or file "
        return
    TargetrMLST = jsonUpGoer(targetjson, markers, target, outdir, 'target')
    typing, removed = compareType(TargetrMLST, nonTargetrMLST)
    json.dump(typing, open(outdir + 'typing.json', 'w'), sort_keys=False, indent=4, separators=(',', ': '))
    json.dump(removed, open(outdir + 'removed.json', 'w'), sort_keys=False, indent=4, separators=(',', ': '))
    for sigtarget in typing:
        SigSeekr(typing, typing[sigtarget], outdir, evalue, float(estop), 200, 1)


#Parser for arguments
parser = ArgumentParser(description='Find Universal Strain-Specifc Probes')
parser.add_argument('--version', action='version', version='%(prog)s v0.5')
parser.add_argument('-o', '--output', required=True, help='Specify output directory')
parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
parser.add_argument('-m', '--marker', required=True, help='Specify rMLST folder')
parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
parser.add_argument('-e', '--evalue', default=1e-1, help='Specify elimination E-value lower limit (default 1e-50)')
parser.add_argument('-s', '--estop', default=1e-70, help='Specify the upper E-value limit (default 1e-90)')
# parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
args = vars(parser.parse_args())

sorter(os.path.join(os.path.abspath(args['marker']), ""),
       os.path.join(os.path.abspath(args['input']), ""),
       os.path.join(os.path.abspath(args['output'])),
       os.path.abspath(args['target']), args['evalue'],
       args['estop'])