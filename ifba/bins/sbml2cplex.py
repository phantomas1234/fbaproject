#!/usr/bin/env python
# encoding: utf-8
"""
sbml2cplex.py

Created by Nikolaus Sonnenschein on 2010-02-22.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
from ifba.GlpkWrap.util import WriteCplex
from pyMetabolism.io.sbml import SBMLParser
from pyMetabolism.fba.glpk import Metabolism2glpk
from pyMetabolism.metabolism import Metabolism

usage = "sbml2cplex.py path2sbml path2cplex"

if sys.argv < 2:
    print usage

parser = SBMLParser(sys.argv[1])
metbol = parser.get_metabolic_system()
converter = Metabolism2glpk(metbol)
lp = converter.convert_to_ifba_metabolism()
WriteCplex(lp, sys.argv[2])


