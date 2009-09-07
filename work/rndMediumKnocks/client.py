#!/usr/bin/env python
# encoding: utf-8
"""
client.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import copy
import ifba.distributedFBA.networking as nw
from ifba.GlpkWrap import util, metabolism, randomMedia, knockouts
import rndMedKnocks as rn

class delegateStub(object):
    def run(self):
        return """Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do 
eiusmod tempor incididunt ut labore et dolore magna
aliqua. Ut enim ad minim veniam, quis nostrud exercitation
ullamco laboris nisi ut aliquip ex ea commodo consequat.
Duis aute irure dolor in reprehenderit in voluptate velit esse
cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat
cupidatat non proident, sunt in culpa qui officia deserunt mollit
anim id est laborum."""


if __name__ == '__main__':
    
    def mainLoop(lp):
        counter = 0
        while True:
            counter += 1
            print counter
            # task = rn.rndMedKnocks(lp, threshold=float(sys.argv[2]))
            task = rn.rndMedKnocks(lp)            
            print "Hello I am the client. I print small caps."
            client = nw.Client(task=task, host=sys.argv[1])
            print "And now my best friend the server."
            client.run()

    # path = '../models/iAF1260template.lp'
    path = '../../ifba/models/iJR904template.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path))
    lp.setReactionObjectiveMinimizeRest('R("R_BiomassEcoli")')
    print lp.getObjective()

    mainLoop(copy.copy(lp))