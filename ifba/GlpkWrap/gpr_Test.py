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
        pass
    
    def testTemplate(self):
        """..."""
        pass

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(test_gpr)
    unittest.TextTestRunner(verbosity=6).run(suite)
