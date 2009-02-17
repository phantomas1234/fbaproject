#!/usr/bin/env python
# encoding: utf-8
"""
singleGeneKnockout.py

Created by Nikolaus Sonnenschein on 2007-11-18.
Copyright (c) 2007 Jacobs University of Bremen. All rights reserved.
"""

import os
import sys
import string
import gzip
import thread
import time
import threading
import Queue
import random
import glpsolWrapper as fba
import util
import string
import re

class Worker(threading.Thread):

    def __init__(self, queue, model):
        self.__queue = queue
        self.model = model
        threading.Thread.__init__(self)

    def run(self):
        print "Hey"
        while 1:
            item = self.__queue.get()
            print item
            if item is None:
                break # reached end of queue
            newObj = self.model.knockOut([item])
            parseObj = fba.GlpsolParse(newObj.glpsol())
            print parseObj
            fluxDist = parseObj.GetFluxDist()
            # if fluxDist[0][1] == 0.:
            #     return (1, fluxDist)
            # else:
            #     return (0, fluxDist)
            if fluxDist[0][1] == 0.:
                print self.getName() ,"  ", item, "finished", fluxDist[0]
            else:
                print self.getName() ,"  ", item, "finished", fluxDist[0]


def main2():
    WORKERS = 4
    virgin = fba.Glpsol(open('test_data/iAF1260template.lp').read())
    count = 1
    reactions = []
    while reactions == []:
        rndMed = virgin.randomMedium()
        parseObj = fba.GlpsolParse(rndMed.glpsol())
        print parseObj
        reactions = parseObj.GetActiveReactions()
    fluxDist = parseObj.GetFluxDist()
    for key in fluxDist[0:10]:
        print 'Reaction:', key[0], 'FluxValue:', key[1]
    print '... viable medium found'
    
    print "starting the singleKnockOut analysis with queue and thread"
    
    targets = util.remTranspRepRev(parseObj.GetActiveReactions())
    print targets
    
    queue = Queue.Queue(20)
    
    for i in range(WORKERS):
        Worker(queue, rndMed).start() # start a worker
    
    for item in targets[1:5]:
        print "push ", item
        queue.put(item)
    
    for i in range(WORKERS):
        queue.put(None) # add end-of-queue markers

if __name__ == '__main__':
    main2()

    