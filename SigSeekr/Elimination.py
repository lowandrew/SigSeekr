#!/usr/bin/env python
from GeneSeekr.GeneSeekrOOP import GeneSeekr
__author__ = 'mike knowles'


def helper(genes, targets, out, cuttoff, threads):
    from glob import glob
    from json import dump
    import time
    import os
    assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
    assert os.path.isfile(genes), u'rMLST fasta not valid {0!r:s}'.format(genes)
    assert isinstance(threads, int)
    ispath = (lambda x: glob(x + "/*.fa*") if os.path.isdir(x) else [x])
    genes = ispath(genes)
    targets = ispath(targets)
    result = GeneSeekr(genes, targets, threads)
    result.mpblast(cuttoff)
    dump(result.plus, open("%s/MLST-gene_results_%s.json" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'w'),
         sort_keys=True, indent=4, separators=(',', ': '))
    result.csvwriter(out, 'rMLST')


class rMLST:

    def _reducer(self, genes):
        overlap = [x for x in self.targets if x not in self.nontargets]
        if len(overlap) is 0:
            rmlst = dict((target, self.nonTargetrMLST[target]) for target in self.targets if target not in overlap)
        if len(overlap) < len(self.targets):
            rmlst = dict((target, self.nonTargetrMLST[target]) for target in self.targets if target not in overlap)
            rmlst.update(GeneSeekr(genes, overlap, self.threads).mpblast(100))
        else:
            rmlst = GeneSeekr(genes, overlap, self.threads).mpblast(100)
        return rmlst

    def __init__(self, targets, nontargets, genes, threads=12):
        self.targets, self.nontargets, self.threads = targets, nontargets, threads
        self.nonTargetrMLST = GeneSeekr(genes, self.nontargets, self.threads).mpblast(100)
        self.TargetrMLST = self._reducer(genes)

    def compare(self):
        typing, removed = dict, dict
        for target in self.targets:
            typing[target], removed[target] = [], []
            for nontarget in self.nontargets:
                similarity = [target, cmp(self.TargetrMLST[target], self.nonTargetrMLST[nontarget])]
                if similarity[1] == 0:
                    removed[target].append([target, similarity])
                else:
                    typing[target].append([target, abs(similarity)])
        return typing, removed