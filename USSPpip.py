__author__ = 'mikeknowles'
"""
Strain-specific probe idenification through:
BLASTing at different e-values
Masking of the resultant
"""
import time, Queue, threading, cStringIO, os, re, sys
from GeneSeekr import makedbthreads
from shutil import copy
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO, SeqIO
from copy import deepcopy
from math import log10
from termcolor import colored

blastqueue = Queue.Queue()
threadlock = threading.Lock()
recorddict = {}
wordsize = [30, 25, 20, 20, 20, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 4, 4, 4, 4]
result = True
minLength = 200


def blastparse(stdout, output, tname, ntname):
    # TODO: ADD minLength variable
    global recorddict, minLength
    handle = open(output, 'w')  # open the target fasta file for writing
    blast_handle = cStringIO.StringIO(stdout)  # Convert string to IO object for use in SearchIO using StringIO
    try:  # Necessary to avoid bad genomes
        for qresult in SearchIO.parse(blast_handle, 'blast-tab'):  # Parse the blast output sting as if it were a file
            for hit in qresult:  # Hit object
                for hsp in hit:  # Hsp object
                    begin = hsp.query_range[0]  # Start of hsp
                    finish = hsp.query_range[1]  # End of hsp
                    if hsp.query_id in recorddict:
                        # For the Contig name in the target fasta dictionary mask using coordinates
                        if finish > begin:
                            recorddict[hsp.query_id].seq = \
                                recorddict[hsp.query_id].seq[:begin] + 'N' * (finish - begin + 1) \
                                + recorddict[hsp.query_id].seq[finish:]
                        else:
                            recorddict[hsp.query_id].seq \
                                = recorddict[hsp.query_id].seq[:finish] + 'N' * (begin - finish + 1) \
                                + recorddict[hsp.query_id].seq[begin:]
        recorddict_bak = deepcopy(recorddict)  # Copy the dictionary so we may iterate and modify the result
        for idline in recorddict_bak:
            # pattern = r'[^N]{'+ re.escape(str(minLength))+r'}' #  Find a sequence of at least the target length
            pattern = r'[ATCG]{30,}N{200,900}[ATCG]{30,}|[^N]{' + re.escape(str(minLength))+r'}'
            if re.match(pattern, str(recorddict[idline].seq)) is not None:
                SeqIO.write(recorddict[idline], handle, "fasta")
            else:
                # print 'Contig \'%s\' not written to file' % id
                recorddict.pop(idline)
    except ValueError:
        print 'Value Error: There was an error removing %s genome from %s' % (ntname, tname)
import pygame
pygame.init()
from pygame.mixer import music


class RunBlast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)

    def run(self):
        while True:
            target, nontarget, tname, ntname, evalue, word, t, count = self.blastqueue.get()
            threadlock.acquire()
            if os.path.getsize(target) != 0:
                # BLASTn parameters
                blastn = NcbiblastnCommandline(query=target,
                                               db=nontarget,
                                               evalue=evalue,
                                               outfmt=6,
                                               perc_identity=90,
                                               word_size=11,
                                               ungapped=True,
                                               penalty=-6)
                stdout, stderr = blastn()
                sys.stdout.write('[%s] [%s/%s] ' % (time.strftime("%H:%M:%S"), count, t))
                if not stdout:
                    print colored("%s has no significant similarity to \"%s\" with an elimination E-value: %g"
                                  % (tname, ntname, evalue), 'red')

                else:
                    music.load('/run/media/blais/blastdrive/coin.wav')
                    music.play()
                    print colored("Now eliminating %s sequences with significant similarity to \"%s\" with an "
                                  "elimination E-value: %g" % (tname, ntname, evalue), 'blue', attrs=['blink'])
                    blastparse(stdout, target, tname, ntname)
            else:
                global result
                result = None
            threadlock.release()
            self.blastqueue.task_done()


def blastnthreads(target, targetnameext, nontargets, targetdir, evalue, word):
    """Setup and create  threads for blastn and passthrough"""
    targetpath = targetdir + targetnameext
    t = len(nontargets)
    count = 0
    if os.path.getsize(targetpath) != 0:
        for nontarget in sorted(nontargets, key=len, reverse=True):
            ntname = nontarget.split('/')[-1].split('.')[0].replace('_', ' ')
            count += 1
            nontargetdb = nontarget.split('.')[0]
            if target != ntname:
                blastqueue.put((targetpath, nontargetdb, target, ntname, evalue, word, t, count))
    else:
        return None
    blastqueue.join()


def restart(target, unique):
    """Write the target fasta file to memory as a dictionary and replace the old target file for iterative purposes"""
    global recorddict
    recorddict = {}
    copy(target, unique)
    # shutil.copy(targetdir + target + '.fa', unique)
    handle = open(target)
    recorddict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    handle.close()


def SigSeekr(targets, nontargets, nontargetdir, evalue, estop, minlength, iterations):
    """Targets and nontargets are imported as a list"""
    global minLength
    minLength = minlength
    # nontargets = []
    unique = nontargetdir + "../Unique/"
    if not os.path.exists(unique):
        os.mkdir(unique)
    makedbthreads(nontargets)
    print "[%s] There are %s target genomes and %s non-target genomes" \
          % (time.strftime("%H:%M:%S"), len(targets), len(nontargets))
    counter = 0
    for target in targets:
        restart(target, unique)
        targetnameext = target.split('/')[-1]
        targetname = targetnameext.split('.')[0]
        counter += 1
        print "[%s] Now parsing target #%i: %s" % (time.strftime("%H:%M:%S"), counter, targetname)
        for i in range(len(nontargets)):
            threads = RunBlast(blastqueue)
            threads.setDaemon(True)
            threads.start()
        for inc in range(iterations):
            print "Iteration " + str(inc + 1)
            if inc > len(wordsize):
                inc = -1
            word = wordsize[inc]
            for e in range(int(log10(evalue)), int(log10(estop)), -1):
                global result
                result = True
                evalue = 10 ** e  # Increment evalue
                blastnthreads(targetname, targetnameext, nontargets, unique, evalue, word)  # BLASTn
                sys.stdout.write('[%s] ' % (time.strftime("%H:%M:%S")))
                if result is not None:
                    print 'Found Sequence(s) at E-value: ' + str(evalue)
                    break
                else:
                    print 'Query file is empty'
                    music.load('/run/media/blais/blastdrive/1up.wav')
                    music.play()
                    restart(target, unique)
                    inc = 0
        uniquename = unique + targetname + '.unique.fasta'
        uniquecount = 0
        targetdict = SeqIO.to_dict(SeqIO.parse(target, 'fasta'))
        handle = open(uniquename, 'w')
        for idline in recorddict:
            pattern = r'([^N]{' + re.escape(str(minLength)) + r',})|([ATCG]{30,}N{' \
                      + re.escape(str(minLength)) + r',900}[ATCG]{30,})'
             #  Find a sequence of at least the target length
            regex = re.compile(pattern)
            uniseq = regex.finditer(recorddict[idline].seq.tostring())
            for coor in uniseq:
                uniquecount += 1
                sequence = targetdict[idline].seq[coor.start():coor.end()].tostring()
                handle.write('>usid%04i_%s_%s\n' % (uniquecount, targetname, idline))
                handle.write(sequence + '\n')
        music.load('/run/media/blais/blastdrive/death.wav')
        print 'Writing %i sequence(s) to file' % uniquecount
        music.play()
        time.sleep(4)
        handle.close()