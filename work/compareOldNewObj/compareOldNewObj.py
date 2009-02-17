#!/usr/bin/env python
# encoding: utf-8
"""
compareOldNewObj.py

Created by Nikolaus Sonnenschein on 2008-10-16.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import metabolism, fluxdist, randomMedia, util
import copy
import matplotlib.pyplot as pyplot


def fluxDistNewOldObj():
    path = '../models/iAF1260template2.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path))
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
    rndlp = randomMedia.Almaas(copy.copy(lp), alwaysInc=include)
    f = rndlp.generateFluxdist()
    print f.getActiveFluxDist()
    condDict = rndlp.currDict
    print f.getFluxDict()['R("R_Ec_biomass_iAF1260_core_59p81M")']
    lp.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    lp.modifyColumnBounds(condDict)
    f2 = lp.fba()
    print f2.getFluxDict()['R("R_Ec_biomass_iAF1260_core_59p81M")']
    return (f.getFluxDict(), f2.getFluxDict())

def scatterData():
    """docstring for scatterData"""
    f, f2 = fluxDistNewOldObj()
    flist = list()
    flist2 = list()
    for key, value in f.iteritems():
        flist.append(value)
        flist2.append(f2[key])
    return (flist, flist2)
    
def makePlot(howMany=100):
    """docstring for plot"""
    for i in range(howMany):
        print i
        data = scatterData()
        fig1 = pyplot.figure(1)
        pyplot.scatter(data[0], data[1])
        pyplot.xlabel('old objective')
        pyplot.ylabel('new objective')
        pyplot.savefig('/Users/niko/tmp/'+'fluxScatter'+ str(i) +'.png')

makePlot()
    