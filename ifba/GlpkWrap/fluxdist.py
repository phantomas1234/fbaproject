#!/usr/bin/env python
# encoding: utf-8
"""
fluxdist.py

Created by Nikolaus Sonnenschein on 2008-02-18.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import gzip
import time
import array
import numpy

def gzipDump(outPutStr, head, ending):
    """docstring for gzipDump"""
    path = head + "%f" % time.time() + ending + '.gz'
    gzip.open(path, 'w').write(outPutStr)

class FluxDist(object):
    """docstring for FluxDist"""
    def __init__(self, lp, prec=1.e-8):
        self.reactions = lp.getColumnIDs()
        self.prec = prec
        self.fluxes = [self._chop(i) for i in lp.primalValues()]
        self.lp = lp

    def _chop(self, value):
        if -1 * self.prec < value < self.prec:
            return 0
        else:
            return value
        
    def __str__(self):
        return str(self.getFluxDict())
        
    def __getitem__(self, key):
        return self.getFluxDict()[key]
        
    def getFluxDict(self):
        return dict(zip(self.reactions, self.fluxes))

    def getActiveFluxDist(self):
        return [i for i in zip(self.reactions, self.fluxes) 
                if self._chop(i[1]) > 0. or self._chop(i[1]) < 0.]
        # CHANGED: Needs a better threshold

    def getActiveReactions(self):
        return [i[0] for i in self.getActiveFluxDist()]
        
    def getActiveColumns(self):
        return [i[0] for i in zip(range(1, len(self.fluxes) + 1), self.fluxes) 
                if i[1] > 0. or i[1] < 0.]
    
    def getActiveMetabolites(self):
        """docstring for getActiveMetabolites"""
        acMetabolites = list()
        acCols = self.getActiveColumns()
        for i in acCols:
            colCoeff = self.lp.getColumnCoef(i)
            for index, value in enumerate(colCoeff):
                if (value != 0.0):
                    metabolite = self.lp.translateRowIndices([index + 1])[0]
                    acMetabolites.append(metabolite)
        return list(set(acMetabolites))
        
    def getFluxArray(self):
        """Returns a array of all fluxes"""
        return numpy.array(self.fluxes, dtype=numpy.float32)
        
    def tsv(self):
        string = str()
        for i in zip(self.reactions, self.fluxes):
            string = string + "%s\t%.16f\n" % i
        return string

    def active_tsv(self):
        string = str()
        for i in self.getActiveFluxDist():
            string = string + "%s\t%.16f\n" % i
        return string


if __name__ == '__main__':
    from ifba.GlpkWrap import util, metabolism
    def main():
        lp = metabolism.Metabolism(util.ImportCplex("./test_data/model.lp"))
        lp.simplex()
        f = FluxDist(lp)
        # print f.getFluxDict()['R("R_BiomassEcoli")']
        # print f.getActiveFluxDist()
        # print f.getFluxArray()
        # print f.getActiveReactions()
        # print f.getActiveColumns()
        print f.getActiveMetabolites()
        print f.active_tsv()

    main()
