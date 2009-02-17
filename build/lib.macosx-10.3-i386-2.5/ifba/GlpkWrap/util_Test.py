#!/usr/bin/env python
# encoding: utf-8
"""
test_util.py

Created by Nikolaus Sonnenschein on 2008-02-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import unittest
import util
from ifba.glpki import glpki
from ifba.GlpkWrap.metabolism import Metabolism


class testImporter(unittest.TestCase):
    def setUp(self):
        pass

class testCplexImporter(testImporter):
    def setUp(self):
        # self.lp = util.ImportCplex('../models/iJR904template.lp')
        # print len(Metabolism(self.lp).getReactions())
        self.lp = util.ImportCplex('./test_data/model.lp')
    
    def testNumCols(self):
        self.assertEqual(glpki.glp_get_num_cols(self.lp), 1473)
    
    def testNumRows(self):
        self.assertEqual(glpki.glp_get_num_rows(self.lp), 904)

# class testFreeMPSImporter(testImporter):
#   def setUp(self):
#       self.lp = util.ImportFreeMPS('test_data/model.mps')
# 
#   def testNumCols(self):
#       self.assertEqual(glpki.lpx_get_num_cols(self.lp), 1606)
# 
#   def testNumRows(self):
#       self.assertEqual(glpki.lpx_get_num_rows(self.lp), 904)


if __name__ == '__main__':
    unittest.main()