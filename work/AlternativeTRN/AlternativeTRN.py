#!/usr/bin/env python
# encoding: utf-8
"""
AlternativeTRN.py

Created by Nikolaus Sonnenschein on 2009-11-12.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import copy
from ifba.GlpkWrap import util, metabolism, randomMedia, fluxdist

def main():
    pass


if __name__ == '__main__':
    path = '../../ifba/models/iAF1260templateMinMax.lp'
    transp = ['R("Mo2b_Transp")', 'R("MglcDb_Transp")', 'R("MlacDb_Transp")']
    lp = metabolism.Metabolism(util.ImportCplex(path))
    # Remove Glucose and Oxygen from the medium
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,0), 'R("MglcDb_Transp")':(0,0)})
    lp.eraseHistory()
    # Test glucose conditions
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,10), 'R("MglcDb_Transp")':(0,10)})
    lp.simplex()
    print lp.getObjVal()
    lp.initialize()
    # Test lactose conditions
    print [(elem, lp.getColumnBounds()[elem]) for elem in transp]
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,10), 'R("MlacDb_Transp")':(0,10)})
    lp.simplex()
    print lp.getObjVal()