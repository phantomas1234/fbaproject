#!/usr/bin/env python
# encoding: utf-8
"""
parallelFBA.py

Created by Nikolaus Sonnenschein on 2009-07-05.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
import pp
from ifba.GlpkWrap import util, metabolism
# from ifba.glpki.glpki import *


def power(a):
    return a**2

def smallPPtest():
    s = pp.Server(ppservers=("212.201.48.108:50000", ), secret="hundekot")
    task = s.submit(power,(99999999,))
    return task


def main():
    pass


if __name__ == '__main__':
    # print power(100)
    # task = smallPPtest()
    # print task()

    lp = metabolism.Metabolism(util.ImportCplex('../models/iAF1260template_minimalMed_noCarb.lp'))
    lp.simplex()
    print lp.getObjVal()
    # s = pp.Server(ncpus=0, ppservers=("212.201.48.108:50000", ), secret="hundekot")
    s = pp.Server()
    job1 = s.submit(lp.getReactions, globals=globals())
    print job1()
    # job2 = s.submit(lp.simplex, depfuncs=(ifba.GlpkWrap.fluxdist.FluxDist, ifba.GlpkWrap.glpk.glp_simplex), modules=('ifba.GlpkWrap.glpk', 'ifba.GlpkWrap.fluxdist'), globals=globals())
    print globals()
    job2 = s.submit(lp.fba)
    print job2()
