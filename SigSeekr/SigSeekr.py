#!/usr/bin/env python
from math import log10
from multiprocessing import Pool

import colored as colored
import sys
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from GeneSeekr.ARMI_Lt import MakeBlastDB, _pickle_method, _unpickle_method
import os
import time

__author__ = 'mike knowles'

class GenObject(object):
    """Object to store static variables"""
    def __init__(self):
        super(GenObject, self).__setattr__('datastore', {})

    def __getattr__(self, key):
        return self.datastore[key]

    def __setattr__(self, key, value):
        self.datastore[key] = value


class Subtractor:
    @staticmethod
    def makeblastdb(target):
        if not os.path.isfile('{}.nhr'.format(target.path)):  # add nhr for searching
            assert os.path.isfile(target.full)  # check that the fasta has been specified properly
            MakeBlastDB('/usr/local/bin/makeblastdb', db=fasta, out=db, dbtype='nucl')()  # Use MakeBlastDB above

    def restart(self, target):
        """Write the target fasta file as a dictionary and replace the old target file for iterative purposes
        :param unique: str contain output unique file
        :param target: str of file location
        """
        from shutil import copy
        copy(target, target.unique)
        handle = open(target)
        self.recorddict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        for record in self.recorddict:  # make sequence mutable (is this still neccesary?)
            self.recorddict[record].seq = self.recorddict[record].seq.tomutable()
        handle.close()

    def compare(self, target):
        full = target
        target = GenObject()
        target.full = full
        target.path, target.ext = os.path.splitext(target)
        target.pwd, target.name = os.path.split(target.path)
        target.unique = "{0}/Unique/{1}.unique.{2}".format(target.pwd, target.path, target.ext)
        return target

    def __init__(self, targets, nontargets, estart=1, estop=1e-70, length=200, i=1, threads=12, out=os.getcwd()):
        self.targets, self.threads = targets, threads
        self.estart, self.estop, self.length, self.i = estart, estop, length, i
        self.nontargets = map(self.compare, nontargets)  # remove the file extension for easier globing
        self.upath = self.db = map((lambda x: os.path.split(x)[0] + "/Unique/"), self.targets)
        print "[{}] Creating necessary databases for BLAST".format(time.strftime("%H:%M:%S"))
        Pool(self.threads).map(self.makeblastdb, self.nontargets)
        self.recorddict = {}
        print "[%s] There are %s target genomes and %s non-target genomes" \
          % (time.strftime("%H:%M:%S"), len(self.targets), len(self.nontargets))



    def _subtractor(self, target):
        pass

    def _aligner(self, target, nontarget, evalue):
        blastn = NcbiblastnCommandline(query=target,
                                       db=nontarget,
                                       evalue=evalue,
                                       outfmt=6,
                                       perc_identity=100,
                                       # word_size=word,
                                       ungapped=True,
                                       penalty=-20,
                                       reward=1,
                                       )
        stdout, stderr = blastn()
        sys.stdout.write('[%s] [%s/%s] ' % (time.strftime("%H:%M:%S"), self.counter, t))
        if not stdout:
            print colored("%s has no significant similarity to \"%s\"(%i) with an elimination E-value: %g"
                          % (tname, ntname, allele, evalue), 'red')
        else:
            # music.load('/run/media/blais/blastdrive/coin.wav')
            # music.play()
            print colored("Now eliminating %s sequences with significant similarity to \"%s\"(%i) with an "
                          "elimination E-value: %g" % (tname, ntname, allele, evalue), 'blue', attrs=['blink'])
            blastparse(stdout, target, tname, ntname)
            # evaluehit += 1


    def unique(self):
        self.counter = 0
        for target in self.targets:
            self.counter += 1
            uniquecount = 0
            target = self.compare(target)
            print "[%s] Now parsing target #%i: %s" % (time.strftime("%H:%M:%S"), counter, target.name)
            for nontarget in self.nontargets:
                for iteration in range(self.i):
                    print "Iteration " + str(iteration + 1)
                    for e in range(self.estart, int(log10(self.estop)), -1):
                        while True: # while self.recordict
                            self._aligner(target, nontarget, 10**e)


