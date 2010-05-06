#!/usr/bin/env python
# encoding: utf-8
"""
rpoZmilp.py

Created by Nikolaus Sonnenschein on 2010-03-02.
Copyright (c) 2010 . All rights reserved.
"""

import sys
import os


if __name__ == '__main__':
    from ifba.GlpkWrap.metabolism import Metabolism
    from ifba.GlpkWrap.util import ImportCplex, WriteCplex
    from ifba.fluxVariability.fluxVariability import Variability
    import random
    import pickle
    
    lp = Metabolism(ImportCplex('../../ifba/models/iAF1260template_minimalMed_noCarb.lp'))
    print len(lp.getColumnIDs())
    lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20), 'R("Mo2b_Transp")':(-20, 20)})
    lp.modifyColumnBounds({'R("R_Ec_biomass_iAF1260_core_59p81M")': (1., 100.)})
    lp.eraseHistory()
    if os.path.exists('blocked.pcl'):
        blocked = pickle.load(open('blocked.pcl'))
    else:
        # blocked = Variability(lp).getBlocked(['R("R_Ec_biomass_iAF1260_core_59p81M")', 'R("R_AGPR")'])
        blocked = Variability(lp).getBlocked()
        sys.exit()
        pickle.dump(blocked, open('blocked.pcl', 'w'))
    print blocked
    print len(blocked)
    lp.deleteReactionsFromStoich(blocked)
    lp.eraseHistory()
    ids = list(set(lp.getColumnIDs()) - set(lp.getTransporters()))
    print ids
    print len(ids)
    lp.addNegativeColumnSwitches(ids)
    # lp.addNegativeColumnSwitches(['R("Mo2b_Transp")'])
    lp.modifyColumnBounds({'R("R_Ec_biomass_iAF1260_core_59p81M")': (1., 100.)})
    binColumns = lp.getColumnsOfType('binary')
    objDict = dict()
    for i in binColumns:
        objDict[i] = 1.
    lp.setObjective(objDict)
    WriteCplex(lp, 'debug.lp')
    lp.toggleVerbosity()
    lp.fba()
    try:
        lp.fba(method='intopt')
    except:
        pass
    lp.mipValues()
    for k, v in lp.mipValues().items():
        if v == 0.:
            print k, v
    # print repr(lp)
