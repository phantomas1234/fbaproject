#!/usr/bin/env python
# encoding: utf-8
"""
milp_cntxtFBA.py

Created by Nikolaus Sonnenschein on 2009-08-29.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
import copy
from ifba.GlpkWrap import util, metabolism
from ifba.expressFBA.ContextFBA import ContextFBA, loadReactionData


def fdsa():
    l, b = helper()
    return l + b

def helper():
    l = 1 + 1
    b = 2 + 2
    c = fdsa()
    return l, b


if __name__ == '__main__':
    model_path = '/Users/niko/arbeit/Data/SBMLmodels/HomoSapiens/humanFBA_ATPObjective_MartinsMedium.lp'
    lp = metabolism.Metabolism(util.ImportCplex(model_path))
    activeReacs1 = set(lp.fba().getActiveReactions())
    print lp.getObjVal()
    lp.deleteReactions(('R("R_GHMT2r")','R("R_GHMT2r_Rev")', 'R("R_GHMT2rm")', 'R("R_GHMT2rm_Rev")'))
    activeReacs2 = set(lp.fba().getActiveReactions())
    print lp.getObjVal()
    print activeReacs1.difference(activeReacs2)

