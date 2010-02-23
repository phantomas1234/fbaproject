#!/usr/bin/env python
# encoding: utf-8
"""
CheckFluxCoupling.py

Created by Nikolaus Sonnenschein on 2009-12-03.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import pickle
from ifba.GlpkWrap import metabolism, util, randomMedia, fluxdist, glpk
from ifba.fluxVariability.fluxVariability import prepareLP, FluxCoupling

usage = "CheckFluxCoupling.py path2model idOfBiomassReaction Reaction1 Reaction2 ..."

if len(sys.argv) < 4:
    print usage
    sys.exit(-1)

printUncoupledFlag = True
lp = prepareLP(metabolism.Metabolism(util.ImportCplex(sys.argv[1], terminal="OFF")), biomRxn=sys.argv[2])
fcLp = FluxCoupling(lp)
fcLp.fluxCouplingFinder(rxns=sys.argv[3:], printUncoupled=printUncoupledFlag)
