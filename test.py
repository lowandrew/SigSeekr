from Bio.Blast.Applications import NcbiblastnCommandline
# from Bio.Blast import NCBIStandalone
from Bio import SearchIO
import cStringIO, time


genome = '/nas/Pipeline_development/USS_Ppip/Genomes/OLC795.fa'

fasta = '/nas/Pipeline_development/rMLST/BACT000001.fas'
blastn = NcbiblastnCommandline(query=genome,
                               db=fasta,
                               outfmt=5,
                               perc_identity=90)
stdout, stderr = blastn()
start = time.time()
blast_handle = cStringIO.StringIO(stdout)

for qresult in SearchIO.parse(blast_handle, 'blast-xml'):
    for hit in qresult:
        for hsp in hit:
            print hsp.query_id, hsp.query_range, hsp.hit_range
            continue
end = time.time() - start
print "Elapsed time %s" % end