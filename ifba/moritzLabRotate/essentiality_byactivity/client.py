#!/usr/bin/env python
# encoding: utf-8
"""
client.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
Modified by Moritz Beber on 2008-11-26.
"""

import sys
import copy
import ifba.distFBA.networking as nw
from ifba.GlpkWrap import util, metabolism, randomMedia, knockouts
import client_task as ct


if __name__ == '__main__':
    
    def mainLoop(lp):
        counter = 0
        while True:
            counter += 1
            print counter
            task = ct.predefMediumFBA(lp)
            print "Hello I am the client. I print small caps."
            client = nw.Client(task=task, host=sys.argv[1])
            print "And now my best friend the server."
            client.run()

    path = '../../models/iAF1260templateMinMax.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path))
    print lp.getObjective()
    mainLoop(lp)

