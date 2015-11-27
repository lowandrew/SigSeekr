__author__ = 'mikeknowles'
from Bio import SeqIO
import re

for record in SeqIO.parse("/nas/Pipeline_development/SPAdesPipelineSandbox/Genome.fasta", "fasta"):
    filename = re.sub('\S*contig\S*|strain|/|\d+_\d+_\d+|\d+_\d+|scaffold\S*|E_coli\S*|genomic.*|\S*-supercont.*'
                      '','', record.description.split('|')[4].split(',')[0].lstrip(), flags=re.IGNORECASE)
    filename = re.sub('_c\d+|\schromosome.*|\sstr|\ssp|\.|\sgenome|\s$', "", filename)
    filename = re.sub('\s$', "", filename)
    filename = re.sub('\||\*|/|\s\s|\s|-|:','_', filename)
    handle = open("/nas/Pipeline_development/SPAdesPipelineSandbox/References/" + filename + '.fasta', 'w')
    handle.write(record.format("fasta"))
    handle.close()
    print filename