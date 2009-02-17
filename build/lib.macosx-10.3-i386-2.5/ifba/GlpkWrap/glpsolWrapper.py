#!/usr/bin/env python
# encoding: utf-8
"""
glpsolWrapper.py

Created by Nikolaus Sonnenschein on 2008-02-08.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import re
import subprocess as subP
import random

class Glpsol:
    """docstring for glpsol"""
    
    def __init__(self, model):
        self.model = model + '\n\n'
    
    def __str__(self):
        """docstring for __str__"""
        print self.model
    
    def getAvailTransp(self, postfix='_Transp'):
        """docstring for getAvailTransp"""
        patt1 = re.compile(r'R\(".+' + postfix + '"\)')
        return list(set(patt1.findall(self.model)))
    
    def getAvailReactions(self):
        """docstring for getAvailReactions"""
        patt2 = re.compile(r'R\("R_.+?"\)')
        return list(set(patt2.findall(self.model)))
    
    def glpsol(self):
        """Pipes a cplex string to the commandline interface of glpk"""
        command = ["glpsol","--cpxlp","/dev/stdin", "-o", "/dev/stdout"]
        (out, err) = subP.Popen(command, stdin=subP.PIPE, stdout=subP.PIPE,
                           stderr=subP.PIPE).communicate(self.model)
        return out
    
    def randomMedium(self, howMany=range(10,100,10), bound=20):
        """docstring for randomMediumFBA"""
        availTransp = self.getAvailTransp()
        l = len(availTransp)/100
        list = [l*i for i in howMany]
        transp = random.sample(availTransp, random.choice(list))
        modelTmp = self.model
        for t in transp:
            lowBound = -1*random.uniform(0, bound)
            highBound =  random.uniform(0, bound)
            bounds = (lowBound, t, highBound)
            modelTmp = modelTmp + '%s <= %s <= %s\n' % bounds
        return Glpsol(modelTmp)
    
    def knockOut(self, reactions):
        """Takes a list of reactions and constrains them all to zero. Also
        reversible counterparts are taken into account. Generates a new glpsol
        object.
        """
        modelTmp = self.model
        for r in reactions:
            modelTmp = modelTmp + '%s <= %s <= %s\n' % (0, r, 0)
            modelTmp = modelTmp + '%s <= %s <= %s\n' % (0, r[:-2] + '_Rev' + r[-2:], 0)
        return Glpsol(modelTmp)
    
    def addBounds(self, medium):
        """Adds dictionary of bounds, e.g.
        {'R("Mo2b_Transp")':(0, 20), 'R("MglcDb_Transp")':(0, 20)},
        to the cplex model and generates a new glpsol object.
        """
        modelTmp = self.model + "\n\ Some model modifications follow\n"
        for key in medium:
            bound = str(medium[key][0]) + " <=  " + key + " <= " + str(medium[key][1])
            modelTmp = modelTmp + bound + "\n"
        return Glpsol(modelTmp)


class GlpsolParse:
    """Parse Glpsol's Output"""
    patt1 = re.compile(r'Column.+?End\sof\soutput', re.S)
    patt2 = re.compile(r'(R\(\".+?\"\)).+?([\d.e-]+)', re.S)
    
    def __init__(self, arg):
        self.arg = GlpsolParse.patt1.search(arg).group()
    
    def __str__(self):
        """docstring for __str__"""
        print self.arg
    
    def GetFluxDist(self):
        """docstring for parse"""
        return [(name, float(val)) for name, val in GlpsolParse.patt2.findall(self.arg)]
    
    def GetActiveFluxDist(self):
        """docstring for GetActiveFluxDist"""
        dist = self.GetFluxDist()
        return filter(lambda x: x[1] != 0 , dist)
    
    def GetActiveReactions(self):
        """  """
        dist = self.GetActiveFluxDist()
        return [elem[0] for elem in dist]

if __name__ == '__main__':
    def testGlpsol():
    	print Glpsol(open('../tests/model.lp', 'r').read()).glpsol()
    testGlpsol()
