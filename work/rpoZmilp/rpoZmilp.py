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
    from ifba.GlpkWrap.util import ImportCplex
    
    lp = Metabolism(ImportCplex('../../ifba/models/iAF1260template.lp'))
    lp.addColumnSwitches(['R("R_PGK")'])
    print repr(lp)
