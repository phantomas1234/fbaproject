#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Nikolaus Sonnenschein on 2008-02-20.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import random
import ifba.combinatorics.combinatorics as comb

class StubAnalysis(object):
    """docstring for StubAnalysis"""
    def __init__(self, viable, deadly):
        super(StubAnalysis, self).__init__()
        self.viable = viable
        self.deadly = deadly
        self.res_list = range(1, 1000 + 1)
        self.immutable_targets = range(1, 1000 + 1)

    def checkIfDeadly(self):
        return random.choice(['DEADLY', 'VIABLE'])
    
    def checkViable(self):
        """docstring for viable"""
        next_res_list = []
        for medium in self.viable:
            print 'Length of res_list', len(self.res_list)
            for reac in self.res_list:
                result = self.checkIfDeadly()
                if result == 'VIABLE':
                    next_res_list.append(reac)
            self.res_list = next_res_list[:]
            next_res_list = []
        print self.res_list
        
    def checkDeadly(self):
        """docstring for viable"""
        next_res_list = []
        for medium in self.deadly:
            print 'Length of res_list', len(self.res_list)
            for reac in self.res_list:
                result = self.checkIfDeadly()
                if result == 'DEADLY':
                    next_res_list.append(reac)
            self.res_list = next_res_list[:]
            next_res_list = []    
        print self.res_list

    def nextLevel(self):
        """docstring for nextLevel"""
        self.res_list = list(comb.SetCombine(self.res_list, self.immutable_targets).generate())

def main():
    o = StubAnalysis([1,2,3,4,5], [1,2,3])
    for i in range(1, 10 + 1): # num of levels
        print "Level: ", i
        o.checkViable()
        o.checkDeadly()
        o.nextLevel()
        
if __name__ == '__main__':
    main()

