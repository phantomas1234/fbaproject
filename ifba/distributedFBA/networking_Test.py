#!/usr/bin/env python
# encoding: utf-8
"""
networking_Test.py

Created by Nikolaus Sonnenschein on 2008-05-8.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import unittest
from ifba.distFBA import networking

class test_Networking(unittest.TestCase):
    def setUp(self):
        pass
        
    def testPriorityQueue(self):
        """Tests the Priority Queue from the Python Cookbook"""
        q = networking.PriorityQueue()
        q.put("2", 2)
        q.put("1", 1)
        q.put("3", 3)
        self.assert_(q.get(), '1')
        self.assert_(q.get(), '2')
        self.assert_(q.get(), '3')

if __name__ == '__main__':
    unittest.main()