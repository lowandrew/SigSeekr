#!/usr/bin/env python
from multiprocessing import Pool
from Bio import SeqIO
from GeneSeekr.ARMI_Lt import MakeBlastDB, _pickle_method, _unpickle_method
import os
import time

__author__ = 'mike knowles'

class GenObject(object):

    def __init__(self):
        super(GenObject, self).__setattr__('datastore', {})

    def __getattr__(self, key):
        return self.datastore[key]

    def __setattr__(self, key, value):
        self.datastore[key] = value


class Subtractor:
    @staticmethod
    def makeblastdb((fasta, db)):
        if not os.path.isfile('{}.nhr'.format(db)):  # add nhr for searching
            assert os.path.isfile(fasta)  # check that the fasta has been specified properly
            MakeBlastDB('/usr/local/bin/makeblastdb', db=fasta, out=db, dbtype='nucl')()  # Use MakeBlastDB above

    def restart(self, target, unique):
        """Write the target fasta file as a dictionary and replace the old target file for iterative purposes
        :param unique: str contain output unique file
        :param target: str of file location
        """
        from shutil import copy
        copy(target, unique)
        handle = open(target)
        self.recorddict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        for record in self.recorddict:  # make sequence mutable (is this still neccesary?)
            self.recorddict[record].seq = self.recorddict[record].seq.tomutable()
        handle.close()

    def __init__(self, targets, nontargets, estart=1, estop=1e-70, length=200, i=1, threads=12, **kwargs):
        self.target = GenObject()
        self.targets, self.nontargets, self.threads = targets, nontargets, threads
        self.db = map((lambda x: os.path.splitext(x)[0]), nontargets)  # remove the file extension for easier globing
        print "[{}] Creating necessary databases for BLAST".format(time.strftime("%H:%M:%S"))
        Pool(self.threads).map(self.makeblastdb, zip(self.nontargets, self.db))
        self.recorddict = {}

    def compare(self, target):
        self.target = GenObject()
        self.target.path, self.target.ext = os.path.splitext(target)
        self.target.pwd, self.target.name = os.path.split(self.target.path)

