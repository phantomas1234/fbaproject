#!/usr/bin/env python
# encoding: utf-8
"""
combinatoricalMedia.py

Created by Nikolaus Sonnenschein on 2008-06-03.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.combinatorics import combinatorics as comb
from ifba.GlpkWrap import metabolism, util, fluxdist
import copy
import itertools
import pickle
import re

def readMinMed():
    f = open('minimalMedium.txt', 'r')
    med = eval(f.read())
    f.close()
    return dict(med)
    
def init(path):
    """A template gets initialized."""
    struct = util.ImportCplex(path)
    blub = metabolism.Metabolism(struct)
    blub.toggleVerbosity()
    return blub
    
def makeMinimalModel():
    med = readMinMed()
    print med
    lp = init('../models/iAF1260template.lp')
    lp.modifyColumnBounds(med)
    return copy.copy(lp)
    
def readSources():
    # f = open('sources.tsv', 'r')
    f = open('sourcesSubSet.tsv', 'r')
    sources = list()
    for line in f:
        line = line.replace('\n', '')
        sources.append(line.split('\t'))
    return [tuple(sub) for sub in sources]
    
def sourceGenerator():
    carbon, nitrogen, phosphate, sulfate =  readSources()
    return ((a, b, c, d) for a in carbon for b in nitrogen for c in phosphate for d in sulfate)
    
def sourceTuple2mediumDict(tup):
    d = dict()
    for i in tup:
        d[i] = (0, 20)
    return d
    
def cleanUpTargets(listOfReactions):
    def revMatch(str):
        patt = re.compile("(_Transp)|(R_EX)")
        if patt.search(str):
            return False
        else:
            return True
    return [r for r in listOfReactions if revMatch(r)]
    
def chop(float):
    if float <= 1e-6 and float >= -1e-6:
        return 0.
    else:
        return float
        

def main():
    lp = makeMinimalModel()
    print lp
    gen = sourceGenerator()
    counter = 0
    mediaGrowth = list()
    # for sourceTuple in itertools.islice(gen, 100):
    for sourceTuple in gen:
        counter += 1
        print counter
        medium = sourceTuple2mediumDict(sourceTuple)
        print sourceTuple
        lp.modifyColumnBounds(medium)
        lp.simplex()
        objValTmp = lp.getObjVal()
        print objValTmp
        if objValTmp > 0.2:
            mediaGrowth.append(sourceTuple)
        lp.initialize()
    print counter
    print len(mediaGrowth)
    f = open('mediaGrowthSmall.txt', 'w')
    pickle.dump(mediaGrowth, f)
    f.close()
    
def main2():
    import pprint
    import array
    f = open('mediaGrowthSmall.txt', 'r')
    sourceTuples = pickle.load(f)
    f.close()

    lp = makeMinimalModel()
    targets = cleanUpTargets(lp.getColumnIDs())
    print len(targets)
    gen = sourceGenerator()
    counter = 0
    
    tmpFile = open('tmpResult.tsv', 'a')
    
    matrix = list()
    # for sourceTuple in itertools.islice(gen, 2):
    for sourceTuple in sourceTuples:
        counter += 1
        print "still ", len(sourceTuples) - counter, " to go."
        medium = sourceTuple2mediumDict(sourceTuple)
        print sourceTuple
        lp.modifyColumnBounds(medium)
        lp.simplex()
        objValTmp = lp.getObjVal()
        vec = list()
        for targ in targets:
            lp.modifyColumnBounds({targ : (0, 0)})
            lp.simplex()
            mutObjVal = chop(lp.getObjVal())
            vec.append(mutObjVal)
            lp.undo()
        lp.initialize()
        matrix.append(vec)
        tmpFile.write('\t'.join([str(i) for i in vec]) + "\n")
        
    tmpFile.close()

    f = open('medKOmatrix.pic', 'wb')
    pickle.dump((matrix, targets, sourceTuple), f)
    f.close()
    
    f = open('medKOmatrix.pic', 'rb')
    struct = pickle.load(f)
    f.close()
    
    pprint.pprint(struct[0])

if __name__ == '__main__':
    main2()

