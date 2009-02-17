#!/usr/bin/env python
# encoding: utf-8
"""
randomMedia_Test.py

Created by Nikolaus Sonnenschein on 2008-02-20.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import unittest
import copy

from ifba.GlpkWrap.metabolism import Metabolism
from ifba.GlpkWrap.randomMedia import Almaas
from ifba.GlpkWrap import fluxdist, util

class randomMedia_Test(unittest.TestCase):
    def setUp(self):
        self.lp = util.ImportCplex('../models/iAF1260template.lp')
        self.glp = Almaas(Metabolism(self.lp))
        # self.glp.toggleVerbosity() # Verbosity is allready toggled 

    def testGenerateFluxDist(self):
        fDist = self.glp.generateFluxdist()
        self.assert_(isinstance(fDist, fluxdist.FluxDist))
        
if __name__ == '__main__':
    unittest.main()