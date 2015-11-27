__author__ = 'mikeknowles'

from USSPpip import SigSeekr
from argparse import ArgumentParser
from glob import glob
from collections import defaultdict

def sorter(genefolder, output, target, evalue, estop):
    genelist = glob(genefolder + "/*.fasta")
    genelist2 = [(gene, 53) for gene in sorted(genelist) if gene not in target]
    # genedict = defaultdict(list)

    # for gene in sorted(genelist):
    #     genelist2.append((gene, 0))
    #     genedict[gene] = genelist
    #     genedict[gene].remove(gene)
    SigSeekr([target], genelist2, output, evalue, estop, 30, 1)



parser = ArgumentParser(description='Find Universal Gene-Specifc Probes')
parser.add_argument('--version', action='version', version='%(prog)s v0.5')
parser.add_argument('-o', '--output', required=True, help='Specify output directory')
parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
parser.add_argument('-e', '--evalue', default=1, help='Specify elimination E-value lower limit (default 1e-90)')
parser.add_argument('-s', '--estop', default=1e-75, help='Specify the upper E-value limit (default 1e-200)')
parser.add_argument('-t', '--target', required=True, help='Specify target genome or folder')
args = vars(parser.parse_args())

import os
# sorter(args['input'], args['output'], args['evalue'], args['estop'])
sorter(os.path.join(os.path.abspath(args['input']), ""),
       os.path.join(os.path.abspath(args['output']), ""),
       os.path.abspath(args['target']), args['evalue'],
       args['estop'])