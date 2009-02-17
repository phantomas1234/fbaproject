#!/usr/bin/env python
# encoding: utf-8
"""
rndMedKnocks.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
Modified by Moritz Beber on 2008-11-26.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism, randomMedia, fluxdist, knockouts
import copy

class predefMediumFBA(object):

    def __init__(self, lp):
        super(predefMediumFBA, self).__init__()
        self.lp = lp

    def run(self, medium):
        """
        Given the specific dictionary 'medium' this function performs a FBA on the
        iAF1260 model under these conditions. Also performs FBA on knock-outs.
        """
        self.lp.modifyColumnBounds(medium)
        try:
            self.lp.glpSimplex()
        except Exception, err:
            print >> sys.stderr, err
            return(dict(), dict())
        fd = fluxdist.FluxDist(self.lp)
        wild_type = self.lp.getObjVal()
        # knock-out each of the active reactions
        rxns = fd.getActiveReactions()
        rdict = dict()
        for rxn in rxns:
            #print rxn, 'is knocked out'
            self.lp.deleteReactions([rxn])
            try:
                self.lp.glpSimplex()
                if self.lp.getObjVal() < 0.0:
                    rdict[rxn] = (1, 1)
                elif self.lp.getObjVal() < (0.05 * wild_type):
                    rdict[rxn] = (1, 1)
                else:
                    rdict[rxn] = (1, 0)
                self.lp.undo()
            except Exception, msg:
                print >> sys.stderr, "Knock-out of reaction", rxn, "has no feasable solution! Assuming 0 growth."
                rdict[rxn] = (1, 1)
                self.lp.undo()
        # knock-out each metabolite mapping of data onto active part needed later
        metbs = fd.getActiveMetabolites()
        mdict = dict()
        for metb in metbs:
            if (metb.endswith('c')):
                #print metb, 'is knocked out.'
                self.lp.deleteMetabolites([metb])
                try:
                    self.lp.glpSimplex()
                    if self.lp.getObjVal() < 0.0:
                        mdict[metb] = (1, 1)
                    elif self.lp.getObjVal() < (0.05 * wild_type):
                        mdict[metb] = (1, 1)
                    else:
                        mdict[metb] = (1, 0)
                    self.lp.undo()
                except Exception, msg:
                    print >> sys.stderr, "Knock-out of metabolite", metb, "has no feasable solution! Assuming 0 growth."
                    mdict[metb] = (1, 1)
                    self.lp.undo()
        # reset modifications for next medium conditions
        self.lp.initialize()
        return (rdict, mdict)
    

if __name__ == '__main__':
    # path = '../../models/iAF1260templateMinMax.lp'
    # lp = metabolism.Metabolism(util.ImportCplex(path))
    # lp.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    # #     # test_without_networking(copy.copy(lp))
    # paths = ['/Volumes/Rudi1/Niko/20081016_relativeEssentiality_0.05/dict'+str(i)+".txt" for i in (38, 39)]
    # print paths
    # conditions = [extractCondition(path) for path in paths]
    # test_without_networking_debug(copy.copy(lp), conditions)
    pass

