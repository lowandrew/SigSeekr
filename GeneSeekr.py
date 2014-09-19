__author__ = 'mikeknowles'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from threading import Thread
from Queue import Queue
from collections import defaultdict
import subprocess, os, glob, time, sys, shlex, re, threading, json, mmap

count = 0


def dotter():
    global count
    # if count <= 80:
    sys.stdout.write('.')
        # count += 1
    # else:
    #     sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
    #     count = 0


def makeblastdb(dqueue):
    while True:  # while daemon
        fastapath = dqueue.get() # grabs fastapath from dqueue
        nhr = "%s.nhr" % (fastapath)
        if not os.path.isfile(str(nhr)):
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)))
            dotter()
        dqueue.task_done() # signals to dqueue job is done
        sys.exit()

# Declare queues, list, and dict
dqueue = Queue()
blastqueue = Queue()
parsequeue = Queue()
testqueue = Queue()
plusqueue = Queue()
plusdict = {}
genedict = defaultdict(list)
blastpath = {}
threadlock = threading.Lock()
# blastexist = {}

def makedbthreads(fastas):
    ''' Setup and create threads for class'''
    for i in range(len(fastas)):
        threads = Thread(target=makeblastdb, args=(dqueue,))
        threads.setDaemon(True)
        threads.start()
    for fasta in fastas:
        dqueue.put(fasta)
    dqueue.join() #wait on the dqueue until everything has been processed

def xmlout (fasta, genome):
    gene = re.search('\/(\w+)\.fas', fasta)
    path = re.search('(.+)\/(.+)\/(.+?)\.fa', genome)
    return path, gene


class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            global blastpath, plusdict
            genome, fasta, blastexist = self.blastqueue.get()
            path, gene = xmlout(fasta, genome)
            out = "%s/tmp/%s.%s.xml" % (path.group(1), path.group(3), gene.group(1))
            threadlock.acquire()
            blastpath[out] = {path.group(3): (gene.group(1),)}
            plusdict[path.group(3)] = {gene.group(1): 0}
            threadlock.release()
            if not os.path.isfile(out):
                dotter()
                blastn = NcbiblastnCommandline(query=genome, db=fasta, evalue=1e-40, out=out, outfmt=5)
                stdout, stderr = blastn()
            if not any(blastpath):
                print out
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
            blastqueue.put((genome, fasta, blastexist))
        blastqueue.join()


class blastparser(threading.Thread): # records, genomes):
    def __init__(self, parsequeue):
        self.parsequeue = parsequeue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            global plusdict, genedict
            xml, genomes, mm, num = self.parsequeue.get()
            records = NCBIXML.parse(mm)
            for record in records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        col = 'N'
                        if hsp.identities == alignment.length:
                            col = alignment.title.split('_')[-1]  # MLST type
                        for genome in genomes:
                            for gene in genomes[genome]:
                                threadlock.acquire()  # precaution
                                if genome not in plusdict:
                                    plusdict[genome] = defaultdict(str)
                                if gene not in plusdict[genome]:
                                    plusdict[genome][gene] = 0
                                plusdict[genome][gene] = col
                                threadlock.release()  # precaution for populate dictionary with GIL
            dotter()
            mm.close()

            self.parsequeue.task_done()

def parsethreader(blastpath, genomes):
    global plusdict
    dotter()
    for i in range(len(genomes)):
        threads = blastparser(parsequeue)
        threads.setDaemon(True)
        threads.start()
    progress = len(blastpath)
    for xml in blastpath:
        handle = open(xml, 'r')
        mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        # time.sleep(0.05) # Previously used to combat open file error
        handle.close()
        parsequeue.put((xml, blastpath[xml], mm, progress))
        parsequeue.join()

def blaster(markers, strains, out, name):
    '''
    The blaster function is the stack manager of the module
    markers are the the target fasta folder that with be db'd and BLAST'd against strains folder
    out is the working directory where the blastxml folder will be placed
    name is the partial title of the csv output
    '''
    global count, genedict, blastpath
    #retrieve markers from input
    fastas = glob.glob(markers + "*.fas")
    #retrieve genomes from input
    genomes = glob.glob(strains + "*.fa")
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    #push markers to threads
    makedbthreads(fastas)
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    if os.path.isfile('%sblastxmldict.json' % strains):
        print "[%s] Loading BLAST data from file" % (time.strftime("%H:%M:%S"))
        blastpath = json.load(open('%sblastxmldict.json' % strains))
    else:
        print "[%s] Now performing BLAST database searches" % (time.strftime("%H:%M:%S"))
        # make blastn threads and retrieve xml file locations
        blastnthreads(fastas, genomes)
        json.dump(blastpath, open('%sblastxmldict.json' % strains, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    print "[%s] Now parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    sys.stdout.write('[%s]' % (time.strftime("%H:%M:%S")))
    parsethreader(blastpath, fastas)
    csvheader = 'Strain'
    row = ""
    rowcount = 0
    for genomerow in plusdict:
        row += "\n" + genomerow
        rowcount += 1
        for generow in plusdict[genomerow]:
            if rowcount <= 1:
                csvheader += ',' + generow
            # for plusrow in plusdict[genomerow][generow]:
            row += ',' + str(plusdict[genomerow][generow])
    with open("%s%s_results_%s.csv" % (out, name, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(csvheader)
        csvfile.write(row)
    return plusdict