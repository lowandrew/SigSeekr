__author__ = 'mikeknowles'

from USSPpip import SigSeekr
from argparse import ArgumentParser
from glob import glob
from multiprocessing import Pool, Lock
from time import sleep

genelist = []
lock = Lock()

def sender(gene):
    SigSeekr([gene, ], zip(genelist, range(len(genelist))), args['output'], args['evalue'], args['estop'], 200, 1)
    # print gene , zip(genelist, range(len(genelist))), args['output'], args['evalue'], args['estop'], 200, 1
    return 1


def sorter(genefolder, output, evalue, estop):
    global genelist
    genelist = glob(genefolder + "*.fasta")


    genelist2 = []
    # genedict = defaultdict(list)
    for gene in sorted(genelist):
        genelist2.append((gene, 0))
    #     genedict[gene] = genelist
    #     genedict[gene].remove(gene)
    zip(genelist, range(len(genelist)))
    # for gene in genelist:

    # p = Pool(len(genelist))
    # p.map(sender, genelist)
    gene = "/nas/knowlesm/Collaborations/SND/combined/CFF.fasta"
    SigSeekr([gene, ], zip(genelist, range(len(genelist))), output, evalue, estop, 200, 1)



parser = ArgumentParser(description='Find Universal Gene-Specifc Probes')
parser.add_argument('--version', action='version', version='%(prog)s v0.5')
parser.add_argument('-o', '--output', required=True, help='Specify output directory')
parser.add_argument('-i', '--input', required=True, help='Specify input genome fasta folder')
parser.add_argument('-e', '--evalue', default=10, help='Specify elimination E-value lower limit (default 1e-90)')
parser.add_argument('-s', '--estop', default=0.1, help='Specify the upper E-value limit (default 1e-200)')
args = vars(parser.parse_args())


sorter(args['input'], args['output'], args['evalue'], args['estop'])