from mmap import  mmap
from shutil import copy
import os
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re

def SigWritter(uniquename, target, uniquecount, targetname, evalue):
    targetdict = SeqIO.to_dict(SeqIO.parse(target, 'fasta'))
    copy(target, uniquename + '.' + str(uniquecount))
    handle = open(uniquename, 'a+')
    if os.path.getsize(uniquename) != 0:
        mm = mmap(handle.fileno(), 0, access=ACCESS_READ)
    else:
        mm = handle
    for idline in recorddict:
        pattern = r'([^N]{' + re.escape(str(minLength)) + r',})|([ATCG]{50,}N{' \
                   + re.escape(str(minLength)) + r',900}[ATCG]{50,})'
         #  Find a sequence of at least the target length
        regex = re.compile(pattern)
        uniseq = regex.finditer(recorddict[idline].seq.tostring())
        for coor in uniseq:
            isunique = True
            uniquecount += 1
            sequence = targetdict[idline].seq[coor.start():coor.end()].tostring()
            handle.seek(0)
            for line in handle:
                if sequence in line:
                    print sequence in line
                    isunique = False
            if isunique is True:
                print 'Found Sequence(s) at E-value: ' + str(evalue)
                handle.write('>usid%04i_%g_%s_%s\n' % (uniquecount, evalue, targetname, idline))
                handle.write(sequence + '\n')
            # else:
            #     global evaluehit
            #     evaluehit = False
    print 'Writing %i sequence(s) to file' % uniquecount
    handle.close()
forward = 'ATCGATGGTGCCTTCGGC'
reverse = 'AAAAGCGGGCAAAACAAAAGG'

# reverse = 'ATCGATGGTGCCTTCGGC'
# forward = 'AAAAGCGGGCAAAACAAAAGG'
reverse = (Seq(reverse)).reverse_complement()
print reverse.tostring()


mishmar = SeqIO.to_dict(SeqIO.parse("/run/media/blais/blastdrive/Salmonella/S-mishmar-haemek-NALR-FFFM_filteredAssembled.fasta", "fasta"))
# reverse = Seq(reverse).reverse_complement().tostring()
for idline in mishmar:
    uniseq = re.compile(forward + r'[ATCG]+' + reverse.tostring()).finditer(mishmar[idline].seq.tostring())
    for coor in uniseq:
        sequence = mishmar[idline].seq[coor.start():coor.end()].tostring()
        print idline
        print coor.start(), coor.end()
        print len(sequence)
        print sequence