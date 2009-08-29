#!/usr/bin/env python
# encoding: utf-8
"""
ExportHumanModels.py

Created by Nikolaus Sonnenschein on 2009-07-13.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism


def writeFullMediumModel(file='humanFullMediumModel.lp'):
    lp = metabolism.Metabolism(util.ImportCplex('../HumanFBA/validHumanFBAmodel.lp'))
    transp = lp.getTransporters()
    print transp
    fullMed = dict()
    for t in transp:
        fullMed[t] = (-1000, 20)
    lp.modifyColumnBounds(fullMed)
    lp.fba()
    print lp
    util.WriteCplex(lp, file)
    

if __name__ == '__main__':
    writeFullMediumModel()

