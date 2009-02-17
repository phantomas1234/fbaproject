#!/usr/bin/env python
# encoding: utf-8
"""
combinatorics_Test.py

Created by Nikolaus Sonnenschein on 2008-02-20.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import unittest
import combinatorics

class combinatorics_Test(unittest.TestCase):
    def setUp(self):
        gen = combinatorics.SetCombine(range(1,20),range(10,30)).generate()
    def testGeneratorType(self):
        """Test if instance of SetCombine is a generator"""

    
if __name__ == '__main__':
    unittest.main()