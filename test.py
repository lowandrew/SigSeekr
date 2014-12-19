# from Bio.Blast.Applications import NcbiblastnCommandline
# from Bio import SearchIO, SeqIO
# import cStringIO, time, glob, re, shutil

cmd = 'blastn -outfmt 5 -query /nas/Pipeline_development/USS_Ppip/test/testing/OLC811.fa -db /nas/Pipeline_development/rMLST/BACT000015 -evalue 1e-40 -perc_identity 100'

genome = '/nas/Pipeline_development/USS_Ppip/test/testing/OLC811.fa'
db = '/nas/Pipeline_development/rMLST/BACT000015'
perc = 100
blastn = NcbiblastnCommandline(query=genome, db=db, evalue=1e-40, outfmt=5, perc_identity=perc)

# genome = '/nas/Pipeline_development/Typing/Genomes/'
# for file in glob.glob(genome + "*'*.f*a"):
#     newname = re.sub('\'', "", file)
#     shutil.move(file, newname)
#     print file
#     print newname

# fasta = '/nas/Pipeline_development/rMLST/BACT000001.fas'
# blastn = NcbiblastnCommandline(query=genome,
#                                db=fasta,
#                                outfmt=6,
#                                perc_identity=90)
# stdout, stderr = blastn()
# start = time.time()
# handle = open(genome, 'r')
# output_handle = open(genome + 's', 'w')
# # mm = mmap.mmap(handle.fileno(), 0)
# recorddict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
# # for id in recorddict:
# #     print id, recorddict[id].seq
#
#
#
# blast_handle = cStringIO.StringIO(stdout)
#
# for qresult in SearchIO.parse(blast_handle, 'blast-tab'):
#     for hit in qresult:
#         for hsp in hit:
#             # print hsp.query_id, hsp.query_range, hsp.hit_range
#             begin = hsp.query_range[0]
#             finish = hsp.query_range[1]
#             # seq = recorddict[hsp.query_id].seq
#             recorddict[hsp.query_id].seq = recorddict[hsp.query_id].seq[:begin] + 'N' * (finish - begin + 1) + recorddict[hsp.query_id].seq[finish:]
#             # print recorddict[hsp.query_id].seq[begin:finish]
#             continue
# for id in recorddict:
#     if not recorddict[id].seq == (len(recorddict[id].seq) + 1) * 'N':
#         SeqIO.write(recorddict[id], output_handle, "fasta")
# end = time.time() - start
# print "Elapsed time %s" % end

import re
p = re.compile("[a-z]")
for m in p.finditer('a1b2c3d4'):
    print m.start(), m.group(), m.end(), 'a1b2c3d4'[m.start():m.end()]

