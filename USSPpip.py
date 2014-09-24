__author__ = 'mikeknowles'
"""
Strain-specific probe idenification through:
BLASTing at different e-values
Masking of the resultant
"""
import time, Queue, threading, GeneSeekr, cStringIO, os, shutil
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO, SeqIO

blastqueue = Queue.Queue()
threadlock = threading.Lock()
recorddict = {}

def blastparse(stdout, output):
    global recorddict
    handle = open(output, 'r+')
    if recorddict == {}:
        recorddict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    blast_handle = cStringIO.StringIO(stdout)
    for qresult in SearchIO.parse(blast_handle, 'blast-tab'):
        for hit in qresult:
            for hsp in hit:
                begin = hsp.query_range[0]
                finish = hsp.query_range[1]
                recorddict[hsp.query_id].seq = recorddict[hsp.query_id].seq[:begin] + \
                                               'N' * (finish - begin + 1) + recorddict[hsp.query_id].seq[finish:]
    for id in recorddict:
        if not recorddict[id].seq == (len(recorddict[id].seq) + 1) * 'N':
            SeqIO.write(recorddict[id], handle, "fasta")
        else:
            recorddict.pop(id)

class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            target, nontarget, tname, ntname = self.blastqueue.get()
            threadlock.acquire()
            blastn = NcbiblastnCommandline(query=target, db=nontarget, evalue=1e-40, outfmt=6, perc_identity=100)
            stdout, stderr = blastn()
            if not stdout:
                print "%s has no significant similarity to %s" % (tname, ntname)
            else:
                print "Now eliminating %s sequences with significant similarity to %s" % (tname, ntname)
                blastparse(stdout, target)
            threadlock.release()
            self.blastqueue.task_done()


def blastnthreads(targets, targetdir, nontargetdir):
    '''Setup and create  threads for blastn and xml path'''
    for target in targets:
        targetPath = targetdir + target + ".fa"
        print "[%s] Now parsing %s" % (time.strftime("%H:%M:%S"), target)
        for i in range(len(target)):
            threads = runblast(blastqueue)
            threads.setDaemon(True)
            threads.start()
        for nontarget in targets[target]:
            nontargetPath = nontargetdir + nontarget + ".fa"
            blastqueue.put((targetPath, nontargetPath, target, nontarget))
        blastqueue.join()


def ssPCR(targets, targetdir, nontargetdir):
    '''
    Targets and nontargets are imported as a list
    '''
    # TODO reassociate path with targets and nontargets and test if files exist
    nontargetPath = []
    unique = targetdir + "../Unique/"
    if not os.path.exists(unique):
        os.mkdir(unique)
    for target in targets:
        shutil.copy(targetdir + target + '.fa', unique)
        for nontarget in targets[target]:
            if nontarget not in nontargetPath:
                nontargetPath.append(nontargetdir + nontarget + '.fa')
    GeneSeekr.makedbthreads(nontargetPath)
    print "[%s] There are %s target genomes and %s non-target genomes" % (time.strftime("%H:%M:%S"), len(targets), len(nontargetPath))
    blastnthreads(targets, unique, nontargetdir)

