#!/usr/bin/env python
# encoding: utf-8
"""
gpr_Test.py

Unittest for the gpr module.

Created by Nikolaus Sonnenschein on 2010-07-29.
Copyright (c) 2010 . All rights reserved.
"""

from ifba.GlpkWrap import gpr
import unittest

class test_gpr(unittest.TestCase):
    def setUp(self):
        

    
    def testSimplex(self):
        """Tests if the real glpk simplex function spits out the correct
        objective value for the iJR904 model under glucose minimal medium
        condition."""
        self.glp.simplex()
        obj = glpk.glp_get_obj_val(self.glp.lp)
        self.assertAlmostEqual(obj, 0.9259122)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(test_glpk)
    unittest.TextTestRunner(verbosity=6).run(suite)
