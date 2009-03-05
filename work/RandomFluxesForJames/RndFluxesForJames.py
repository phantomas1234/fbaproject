#!/usr/bin/env python
# encoding: utf-8
"""
RndFluxesForMCanalysis.py

Created by Nikolaus Sonnenschein on 2009-01-23.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import copy
import gzip
from ifba.GlpkWrap import metabolism, randomMedia, util
from multiprocessing import Pool

def helpIterator(lp, r):
    for i in range(r):
        yield (lp, i)

def mainStub(args): return main(*args)

def main(lp, stub):
    include = ('R("Mhb_Transp")',
             'R("Mna1b_Transp")',
             'R("Mkb_Transp")',
             'R("Mca2b_Transp")',
             'R("Mcu2b_Transp")',
             'R("Mmg2b_Transp")',
             'R("Mzn2b_Transp")',
             'R("Mmobdb_Transp")',
             'R("Mfe2b_Transp")',
             'R("Mfe3b_Transp")',
             'R("Mcobalt2b_Transp")',
             'R("Mmn2b_Transp")',
             'R("Mclb_Transp")')
    lp = randomMedia.Almaas(lp, alwaysInc=include)
    f = lp.generateFluxdist()
    stringDump = randomMedia.dict2tsv(lp.currDict) + "\n" + f.active_tsv()
    path = "/Users/niko/tmp/" + "iAf1260_fluxDist_" + str(stub) + ".tsv.gz"
    gzip.open(path, 'w').write(stringDump)
    del(lp)
    return stub

if __name__ == '__main__':
    pool = Pool(processes=2)
    tmpLP = metabolism.Metabolism(util.ImportCplex("../../ifba/models/iAF1260template2.lp"))
    tmpLP.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    result = pool.map(mainStub, helpIterator(tmpLP, 100000))