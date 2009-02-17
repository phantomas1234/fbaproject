#!/usr/bin/env python
# encoding: utf-8
"""
util.py

Created by Nikolaus Sonnenschein on 2008-02-08.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""


import string
import re
from ifba.glpki import glpki

def ImportCplex(file_path, terminal="OFF"):
    """Returns a lp struct which can be us by"""
    if terminal == "OFF":
        glpki.glp_term_out(glpki.GLP_OFF)
    elif terminal == "ON":
        glpki.glp_term_out(glpki.GLP_ON)
    else:
        raise Exception, 'wrong option specified.'
    prob = glpki._glp_lpx_create_prob()
    glpki.glp_read_lp(prob ,None, file_path)
    return prob
    # return glpki._glp_lpx_read_cpxlp(file_path)

def ImportMPS(file_path):
    """Returns a lp struct which can be us by"""
    glpki.glp_term_out(glpki.GLP_OFF)
    return glpki.glp_read_freemps(file_path)

def WriteCplex(lp, file_path='debug.lp'):
    """Returns a lp struct which can be us by"""
    glpki.glp_write_lp(lp.lp, file_path)


# =====================
# = Utility functions =
# =====================

def revMatch(str):
    patt = re.compile(".+_Transp")
    if patt.search(str):
        return False
    else:
        return True

def transporters(reactionList):
    return [elem for elem in reactionList if not revMatch(elem)]

def remTranspRepRev(reactionList):
    """docstring for removeTranspReplaceRev"""
    rTmp = [string.replace(elem, '_Rev', '') for elem in reactionList]
    return [elem for elem in rTmp if revMatch(elem)]

def dict2mathematica(dict):
    """docstring for fluxDict2mathematica"""
    out = string.replace(str(dict), ':', ' ->')
    out = string.replace(out, "'", "")
    out = string.replace(out, "(", "[")
    out = string.replace(out, ")", "]")
    return out

def gzipDump(outPutStr, head):
    """docstring for gzipDump"""
    path = head + "%f" % time.time() + '.m' + '.gz'
    gzip.open(path, 'w').write(outPutStr)

def rndMed(transpList, bound=20, howMany=100):
    transp = random.sample(transpList,howMany)
    newBnds = {}
    for elem in transp:
        lowBound = -1*random.uniform(0, bound)
        highBound =  random.uniform(0, bound)
        newBnds[elem] = (lowBound, highBound)
    return newBnds

if __name__ == '__main__':
    lp = ImportCplex('../models/iAF1260template.lp')
    print glpki.lpx_get_num_cols(lp)
    lp2 = ImportFreeMPS('test_data/model.mps')
    print glpki.lpx_get_num_cols(lp2)

