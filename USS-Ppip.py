__author__ = 'mikeknowles'
"""
Strain-specific probe idenification through:
BLASTing at different e-values
Masking of the resultant
"""
import time, Queue, threading
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIStandalone

blastqueue = Queue.Queue()

class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            global blastpath, plusdict
            genome, fasta, blastexist = self.blastqueue.get()

            blastn = NcbiblastnCommandline(query=genome, db=fasta, evalue=1e-40, outfmt=5, perc_identity=100)
            stdout, stderr = blastn()
            self.blastqueue.task_done()


def blastnthreads(fastas, genomes):
    '''Setup and create  threads for blastn and xml path'''
    blastexist = {}
    for i in range(len(fastas)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        threads.start()
    for genome in genomes:
        for fasta in fastas:
            blastqueue.put((genome, fasta))
        blastqueue.join()


def ssPCR(targets, nontargets, targetdir, nontargetdir):
    '''
    Targets and nontargets are imported as a list
    '''
    # TODO reassociate path with targets and nontargets and test if files exist
    print "[%s] There are %s target genomes and %s non-target genomes"
    print "[%s] Now parsing %s"

