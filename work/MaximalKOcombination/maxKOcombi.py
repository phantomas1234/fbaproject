#!/usr/bin/env python
# encoding: utf-8
"""
maxKOcombi.py

Created by Nikolaus Sonnenschein on 2008-06-02.
Copyright (c) 2008 Jacob University Bremen. All rights reserved.
"""

import sys
import os
from ifba.fluxVariability.fluxVariability import Variability
from ifba.rpoZanalysis.rpoZanalysis import *


def loadTheModel():
    # loading the targets
    deadly = loadSource('../rpoZanalysis/targetCarbonSources.txt')
    # loading the viable sources
    viable = loadSource('../rpoZanalysis/viableCarbonSources.txt')
    # loading minimal medium conditions
    minMed = loadMinimalMed('../rpoZanalysis/minimalMedium.txt')
    lp = init('../models/iAF1260template.lp')
    lp.modifyColumnBounds(minMed)
    # lp.simplex()
    # print lp.getObjVal()
    return copy.copy(lp), deadly, viable

def getVariablity(lp):
    varGenerator = Variability(lp, noise=1.).variabilityGenerator()
    vResult = []
    for var in varGenerator:
        print var
        vResult.append(var)
    return vResult

def main():
    lp, deadly, viable = loadTheModel()
    print lp, deadly, viable
    lp.modifyColumnBounds({'R("MglcDb_Transp")' : (0, 20)})
    lp = copy.copy(lp)
    
    # now check the variablity of every flux (obj val is fixed)
    
    result = getVariablity(lp)
    print [i[1] for i in result]
    open('VariabilityAnalysis.txt', 'w').write(repr(result))
    


if __name__ == '__main__':
	main()

