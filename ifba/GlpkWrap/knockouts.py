#!/usr/bin/env python
# encoding: utf-8
"""
knockouts.py

Created by Nikolaus Sonnenschein on 2008-02-25.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

from metabolism import Metabolism
from util import ImportCplex
from ifba.glpki.glpki import *


class KnockOut(object):
    """A class putting the necessary functionality for GeneKnockOuts to
    the Metabolism class."""
    def __init__(self, lp):
        # super(KnockOut, self).__init__(lp)
        self.lp = lp
        self.lp.smcp.presolve = GLP_OFF
    
    def knockOut(self, gene):
        """Knocks out a gene."""
        self.lp.modifyColumnBounds({gene: (0., 0.)})
    
    def knockOuts(self, listOfGenes):
        """Knocks out a list of genes."""
        for gene in listOfGenes:
            self.knockOut(gene)
        
    
if __name__ == '__main__':

    def init(path):
        struct = ImportCplex(path)
        return Metabolism(struct)

    def main():
        import util
        ecoli = init('test_data/model.lp')

        ecoli.simplex()
        print ecoli.getObjVal()
        
        KnockOut(ecoli).knockOuts(['R("R_PGK")', 'R("R_PGK_Rev")'])
        ecoli.simplex()
        print ecoli.getObjVal()
        
    main()