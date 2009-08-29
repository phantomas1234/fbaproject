#!/usr/bin/env python
# encoding: utf-8
"""
CheckBlocked.py

Created by Nikolaus Sonnenschein on 2009-08-20.
Copyright (c) 2009 . All rights reserved.
"""

from ifba.GlpkWrap import util, metabolism
from ifba.blockedReactions.blockedReactions import blockedQ, analyseBlockedReactions
import csv

def checkBlocked(model_path, reacs=None):
    lp = metabolism.Metabolism(util.ImportCplex(model_path))
    lp.simplex()
    print lp
    if reacs:
        result = analyseBlockedReactions(lp, reacs)
    else:
        result = analyseBlockedReactions(lp)
    return result
    
def writeblocked(path, iterator):
    fileobj = open(path, 'w')
    csv_writer = csv.writer(fileobj, dialect='excel-tab')
    for row in iterator:
        csv_writer.writerow(row)
    fileobj.close()
    

if __name__ == '__main__':
    outputPath = '/Users/niko/arbeit/Publishing/Arndt_Colab/data/blockedReactions/'
    model_path = '/Users/niko/arbeit/Data/SBMLmodels/HomoSapiens/humanFBA_ATPObjective_MartinsMedium.lp'
    # blocked = checkBlocked(model_path, ('R("R_GHMT2r")','R("R_GHMT2r_Rev")'))
    blocked = checkBlocked(model_path)
    writeblocked(outputPath+'atp_max_blocked.tsv',[(i,) for i in blocked])
    model_path = '/Users/niko/arbeit/Data/SBMLmodels/HomoSapiens/humanFBA_ATPObjective_MartinsMedium_noOxygen.lp'
    blocked = checkBlocked(model_path)
    writeblocked(outputPath+'atp_max_noOxygen_blocked.tsv',[(i,) for i in blocked])
    model_path = '/Users/niko/arbeit/Data/SBMLmodels/HomoSapiens/humanFBAbiomassObjective.lp'
    blocked = checkBlocked(model_path)
    writeblocked(outputPath+'biomass_max_blocked.tsv',[(i,) for i in blocked])
    
