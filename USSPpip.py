__author__ = 'mikeknowles'
"""
Strain-specific probe idenification through:
BLASTing at different e-values
Masking of the resultant
"""
import time, Queue, threading, GeneSeekr, cStringIO, os, shutil, re, sys
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

def blastparse(stdout, output, tname, ntname):
    # TODO: ADD minLength variable
    minLength = 60
    global recorddict
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
                            recorddict[hsp.query_id].seq = recorddict[hsp.query_id].seq[:begin] + \
                                                        'N' * (finish - begin + 1) + recorddict[hsp.query_id].seq[finish:]
                        else:
                            recorddict[hsp.query_id].seq = recorddict[hsp.query_id].seq[:finish] + \
                                                        'N' * (begin - finish + 1) + recorddict[hsp.query_id].seq[begin:]
        recorddict_bak = deepcopy(recorddict)  # Copy the dictionary so we may iterate and modify the result
        for id in recorddict_bak:
            pattern = r'[^N]{'+ re.escape(str(minLength))+r'}' #  Find a sequence of at least the target length
            try:
                if (re.search(pattern, str(recorddict[id].seq))).group(0):
                    SeqIO.write(recorddict[id], handle, "fasta")
            except AttributeError:
                # print 'Contig \'%s\' not written to file' % id
                recorddict.pop(id)
    except ValueError:
        print 'Value Error: There was an error removing %s genome from %s' % (ntname, tname)

class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            target, nontarget, tname, ntname, evalue = self.blastqueue.get()
            threadlock.acquire()
            if os.path.getsize(target) != 0:
                # BLASTn parameters
                blastn = NcbiblastnCommandline(query=target,
                                               db=nontarget,
                                               evalue=evalue,
                                               outfmt=6,
                                               perc_identity=99,
                                               word_size=30)
                stdout, stderr = blastn()
                sys.stdout.write('[%s] ' % (time.strftime("%H:%M:%S")))
                if not stdout:
                    print colored("%s has no significant similarity to %s with an elimination E-value: %s" \
                          % (tname, ntname, str(evalue)), 'red')
                else:
                    print colored("Now eliminating %s sequences with significant similarity to %s with an elimination " \
                          "E-value: %s" % (tname, ntname, str(evalue)), 'green', attrs=['blink'])
                    blastparse(stdout, target, tname, ntname)
            else:
                global result
                result = False
                print 'Query file is empty'
            threadlock.release()
            self.blastqueue.task_done()


def blastnthreads(target, targets, targetdir, nontargetdir, evalue):
    '''Setup and create  threads for blastn and xml path'''
    targetPath = targetdir + target + ".fa"
    for i in range(len(target)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        threads.start()
    if os.path.getsize(targetPath) != 0:
        for nontarget in sorted(targets[target]):
            nontargetPath = nontargetdir + nontarget + ".fa"
            blastqueue.put((targetPath, nontargetPath, target, nontarget, evalue))
    else:
        return False
    blastqueue.join()

def restart(target, targetdir, unique):
    '''
    Write the target fasta file to memory as a dictionary and replace the old target file for iterative purposes
    '''
    global recorddict
    recorddict = {}
    shutil.copy(targetdir + target + '.fa', unique)
    handle = open(targetdir + target + '.fa', unique)
    recorddict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    handle.close()


def ssPCR(targets, targetdir, nontargetdir):
    '''
    Targets and nontargets are imported as a list
    '''
    # TODO reassociate path with targets and nontargets and test if files exist
    minLength = 60
    evalue = 1e-100
    nontargetPath = []
    unique = targetdir + "../Unique/"
    if not os.path.exists(unique):
        os.mkdir(unique)
    for target in sorted(targets):
        restart(target, targetdir, unique)
        for nontarget in targets[target]:
            if nontarget not in nontargetPath:
                nontargetPath.append(nontargetdir + nontarget + '.fa')
    GeneSeekr.makedbthreads(nontargetPath)
    print "[%s] There are %s target genomes and %s non-target genomes" % (time.strftime("%H:%M:%S"), len(targets), len(nontargetPath))
    for target in targets:
        print "[%s] Now parsing %s" % (time.strftime("%H:%M:%S"), target)
        for e in range(int(log10(evalue)), -200, -1):
            evalue = 10 ** e
            blastnthreads(target, targets, unique, nontargetdir, evalue)
            sys.stdout.write('[%s] ' % (time.strftime("%H:%M:%S")))
            global result
            if result:
                # TODO: Add iterable word size
                print 'Found Sequence(s) at E-value: ' + str(evalue)
                break
            else:
                print 'Query file is empty'
                restart(target, targetdir, unique)
        uniquename = unique + target + '.unique.fasta'
        uniquecount = 1
        handle = open(uniquename, 'w')
        for id in recorddict:
            pattern = r'([^N]{' + re.escape(str(minLength)) + r',})' #  Find a sequence of at least the target length
            regex = re.compile(pattern)
            uniseq = regex.findall(str(recorddict[id].seq))
            print 'Writing %s sequence(s) to file' % str(len(uniseq))
            for sequence in uniseq:
                handle.write('>usid%s_%s_%s\n' % (uniquecount.zfill(4), target, id))
                handle.write(sequence[1] + '\n')
                uniquecount += 1
        handle.close()
