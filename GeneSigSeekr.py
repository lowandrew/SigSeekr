__author__ = 'mikeknowles'

from USSPpip import SigSeekr
from argparse import ArgumentParser
from glob import glob
from collections import defaultdict

def sorter(genefolder, output, evalue, estop):
    genelist = glob(genefolder + "/*.fasta")
    genelist2 = []
    # genedict = defaultdict(list)
    for gene in sorted(genelist):
        genelist2.append((gene, 0))
    #     genedict[gene] = genelist
    #     genedict[gene].remove(gene)
    SigSeekr(genelist, genelist, output, evalue, estop, 200, 1)



parser = ArgumentParser(description='Find Universal Gene-Specifc Probes')
parser.add_argument('--version', action='version', version='%(prog)s v0.5')
parser.add_argument('-o', '--output', required=True, help='Specify output directory')
parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
parser.add_argument('-e', '--evalue', default=1e-100, help='Specify elimination E-value lower limit (default 1e-90)')
parser.add_argument('-s', '--estop', default=1e-250, help='Specify the upper E-value limit (default 1e-200)')
args = vars(parser.parse_args())


sorter(args['input'], args['output'], args['evalue'], args['estop'])