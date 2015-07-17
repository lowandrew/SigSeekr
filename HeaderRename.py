__author__ = 'mikeknowles'

from glob import glob
from Bio import SeqIO
from multiprocessing import Pool

def renamer(genome):
    genomename = genome.split("/")[-1].split(".")[0]
    handle = open(genome)
    records = SeqIO.parse(handle, 'fasta')
    outfile = open("/nas/Genomes/ePCR/%s.fasta" % genomename, 'w')
    for record in records:
        if "gi|" in record.id:
            record.id = genomename + "_" + record.id.split("|")[1]
        else:
            record.id = genomename
        SeqIO.write(record, outfile, 'fasta')
    handle.close()
    outfile.close()
    return 0

p = Pool(20)
p.map(renamer, glob("/nas/Genomes/Shigella/*.fasta"))
# for genome in glob("/run/media/blais/blastdrive/EcoliAll/*.f*"):
#     renamer(genome)