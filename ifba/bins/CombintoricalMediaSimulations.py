#!/usr/bin/env python
# encoding: utf-8
"""
CombinatoricalMediaSimulations.py

Simulates all possible minimal media compositions consisting of unique carbon,
nitrogen, sulfate, phosphate sources. If a single compound provides multiple
elemental sources it serves a the solely source for them. # TODO: Reformulate this

KO analyses of genes and reactions are optional.
TODO:
1. Yaml has to be replaced with JSON (since 2.6 JSON is part of the stdlib)


Created by Nikolaus Sonnenschein on 2010-11-03.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
import time
import textwrap
import Queue
import yaml
# try:
#     import json
# except ImportError:
#     import simplejson as json
if sys.argv[1] != 'client':
    from ifba.storage.hdf5storage import SimulationDB, h5Container
# except ImportError:
#     
#         pass
#     else:
#         print "HDF5/pyTables is not available! You can use the script only in client mode!"
#         sys.exit(-1)
from ifba.distributedFBA.networking import Server, Client
from ifba.distributedFBA.concurrency import GeneratorInputClient, h5OutputClient
from ifba.GlpkWrap.util import ImportCplex
from ifba.GlpkWrap.metabolism import Metabolism
from ifba.GlpkWrap.fluxdist import FBAsimulationResult
from RandomMediaSimulations import generateStorageObject
from ifba.glpki.glpki import glp_delete_prob

def generateStorageObject(path, lp):
    """docstring for generateStorageObject"""
    return SimulationDB(h5Container(path, lp))

def readSources(path):
    f = open(path, 'r')
    sources = list()
    for line in f:
        line = line.replace('\n', '')
        sources.append(line.split('\t'))
    return [tuple(sub) for sub in sources]

def generateCombinatoricalSets(*sources):
    tmpString1 = " ".join(["for e%d in sources[%d]" % (i, i) for i in range(len(sources))])
    tmpString2 = "("+ ", ".join(["e%d" % i for i in range(len(sources))]) + ")"
    funcString = "gen = ("+tmpString2+" "+tmpString1+")"
    codeBlock = compile(funcString, "blah", "exec")
    exec(codeBlock, locals())
    return gen

def generateCombinatoricalMedia(uptake, setGenerator):
    for s in setGenerator:
        yield dict([('R("'+elem+'_Transp")', (0, uptake)) for elem in s])

def readSourcesTable(path):
    """Reads a table of the form:
            source  carbon  nitrogen    phosphor    sulfur
            carb1   1       1           0           1"""
    auxFunc = lambda x: (x[0], int(x[1]), int(x[2]), int(x[3]), int(x[4]))
    return [auxFunc(line.rstrip().split('\t')) for line in open(path)]

def possibleSigs(sig):
    return [(elem1, elem2, elem3, elem4) \
            for elem1 in [0, 1][0:abs(sig[0] - 1)+1] \
            for elem2 in [0, 1][0:abs(sig[1] - 1)+1] \
            for elem3 in [0, 1][0:abs(sig[2] - 1)+1] \
            for elem4 in [0, 1][0:abs(sig[3] - 1)+1]][1:]

def _DirtyHack(sourceTable, seedSources):
    sources2sig = dict([(row[0], tuple(row[1:]))for row in sourceTable])
    combiDict =  dict([(k, list()) for k in [(elem1, elem2, elem3, elem4) for elem1 in [0, 1] for elem2 in  [0, 1] for elem3 in [0, 1] for elem4 in [0, 1]]])
    for row in sourceTable:
        combiDict[tuple(row[1:])].append(row[0])
    for i, carb in enumerate(seedSources):
        sig = sources2sig[carb]
        # print carb, sig
        # print "possible other signatures:"
        nextSigs = possibleSigs(sig)
        if nextSigs == []:
            yield [carb]
        for sig2 in nextSigs:
            potentialSources = combiDict[sig2]
            # print potentialSources
            tmp1 = [[carb, s] for s in potentialSources]
            # print '\t', sig2
            combSig1 = tuple([sig[i] + sig2[i] for i in range(4)])
            # print "\tpossible other signatures2:"
            nextSigs = possibleSigs(combSig1)
            if nextSigs == []:
                for elem in tmp1:
                    yield elem
            for sig3 in nextSigs:
                potentialSources = combiDict[sig3]
                # print potentialSources
                # if potentialSources == []:
                #     for elem in tmp1:
                #         yield elem
                tmp2 = [elem + [s] for elem in tmp1 for s in potentialSources]
                # print '\t\t', sig3
                combSig2 = tuple([combSig1[i] + sig3[i] for i in range(4)])
                # print "\t\tpossible other signatures2:"
                nextSigs = possibleSigs(combSig2)
                if nextSigs == []:
                    for elem in tmp2:
                        yield elem
                for sig4 in nextSigs:
                    potentialSources = combiDict[sig4]
                    # print "4", potentialSources
                    # if potentialSources == []:
                    #     for elem in tmp2:
                    #         yield elem
                    # print 3*'\t', sig4
                    # print textwrap.fill(str(tmp2), initial_indent=3*'\t', subsequent_indent=3*'\t')
                    for s in potentialSources:
                        for elem in tmp2:
                            yield elem + [s]

def combinatoricalSources(sourceTable, seedSources):
    mediaGenerator = _DirtyHack(sourceTable, seedSources)
    stuff = set()
    for med in mediaGenerator:
        # print 4*'\t', med
        stuff.add(tuple(sorted(med)))
    # print len(stuff)
    return stuff

def testCombinatoricalSources(sourceTable, carbonSources):
    import numpy
    stuff = combinatoricalSources(sourceTable, carbonSources)
    stuff3 = set()
    for elem in list(stuff):
        tot = numpy.array([sources2sig[s] for s in elem]).sum()
        if tot != 4:
            print elem
            print [sources2sig[s] for s in elem]
        stuff3.add(tot)
    print stuff3



def generateSolveMediumObject(path2model="", medium={}, include={}, objective=None, optimizationRoutine='pFBA', koQ=True, *args, **kwargs):
    return SolveMedium(path2model=path2model, medium=medium, include=include, objective=objective, optimizationRoutine=optimizationRoutine, koQ=koQ, *args, **kwargs)

class SolveMedium(object):
    
    def __init__(self, path2model="", medium={}, include={}, objective=None, optimizationRoutine='pFBA', koQ=True, *args, **kwargs):
        self.koQ = koQ
        self.optimizationRoutine = optimizationRoutine
        self.objective = objective
        self.lp = Metabolism(ImportCplex(path2model))
        self.path2model = path2model
        if objective:
            self.lp.setReactionObjective(self.objective)
        self.preMed = dict([(r, (-1000., 0)) for r in self.lp.getTransporters()])
        self.preMed.update(include)
        self.lp.modifyColumnBounds(self.preMed)
        self.lp.modifyColumnBounds(dict([(r, (0., 1000.)) for r in self.lp.getReactions()]))
        self.lp.modifyColumnBounds(medium)
        self.lp.eraseHistory()

    def run(self, *args, **kwargs):
        """docstring for run"""

        f = getattr(self.lp, self.optimizationRoutine)()
        knockoutEffects = dict()
        wt = f[self.objective]
        if self.koQ and wt > 0.:
            knockoutEffects = self.lp.singleKoAnalysis(f.getActiveReactions())
            for k in knockoutEffects:
                knockoutEffects[k] = knockoutEffects[k] / wt
        self.lp.undo()
        return FBAsimulationResult(f, knockoutEffects, self.lp.getColumnBounds(),
                                                self.lp.getObjectiveFunction(),
                                                time.time(), self.path2model, "Test")

    def __del__(self):
        """docstring for __del__"""
        glp_delete_prob(self.lp) # FIXME this is a dirty hack
        del self


def solveMedium(path2model="", medium={}, include={}, objective=None, optimizationRoutine='pFBA', koQ=True, *args, **kwargs):
    """doc"""
    lp = Metabolism(ImportCplex(path2model))
    if objective:
        lp.setReactionObjective(objective)
    preMed = dict([(r, (-1000., 0)) for r in lp.getTransporters()])
    preMed.update(include)
    lp.modifyColumnBounds(preMed)
    lp.modifyColumnBounds(medium)
    lp.modifyColumnBounds(dict([(r, (0., 1000.)) for r in lp.getReactions()]))
    lp.eraseHistory()
    # print lp.cplex()
    f = lp.pFBA()
    # simulationStorage = generateStorageObject(outputfile, lp)
    knockoutEffects = dict()
    wt = f[objective]
    print wt
    if koQ and wt > 0.:
        knockoutEffects = lp.singleKoAnalysis(f.getActiveReactions())
        for k in knockoutEffects:
            knockoutEffects[k] = knockoutEffects[k] / wt
    lp.initialize()
    # print knockoutEffects
    return FBAsimulationResult(f, knockoutEffects, lp.getColumnBounds(),
                                            lp.getObjectiveFunction(),
                                            time.time(), path2model, "Test")

def basicFunctionality(outputfile, configPath):
    config = yaml.load(open(configPath))
    descr = yaml.dump(config)
    print descr
    config['descr'] = descr
    sourceTable = readSourcesTable(config["sourcesPath"])
    carbonSources = [elem[0] for elem in sourceTable if elem[1] == 1]
    sources2sig = dict([(row[0], tuple(row[1:]))for row in sourceTable])
    combSources = list(combinatoricalSources(sourceTable, carbonSources))
    gen = generateCombinatoricalMedia(config["uptake"], combSources)
    run = 0
    for medium in gen:
        print "Run:", run
        run += 1
        solveMedium(medium=medium,**config)


def client(serverip):
    """docstring for client"""
    counter = 0
    
    client = Client(task=generateSolveMediumObject, host=serverip)
    while True:
        counter = counter + 1
        print counter
        client.run()

def stub(gen, config):
    for elem in gen:
        config["medium"] = elem
        yield config

def server(outputfile='test.h5', configPath='parameters.yaml'):
    """Server"""
    config = yaml.load(open(configPath))
    descr = yaml.dump(config)
    print descr
    config['descr'] = descr
    sourceTable = readSourcesTable(config["sourcesPath"])
    carbonSources = [elem[0] for elem in sourceTable if elem[1] == 1]
    sources2sig = dict([(row[0], tuple(row[1:]))for row in sourceTable])
    combSources = list(combinatoricalSources(sourceTable, carbonSources))
    gen = generateCombinatoricalMedia(config["uptake"], combSources)
    lp = Metabolism(ImportCplex(config["path2model"]))
    simulationStorage = generateStorageObject(outputfile, lp)
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    # gen2 = (config["medium"] = elem for elem in gen)
    gen2 = stub(gen, config)
    t1 = GeneratorInputClient(inputQueue, gen2)
    t1.start()
    time.sleep(1)
    t2 = h5OutputClient(outputQueue, simulationStorage)
    t2.start()
    time.sleep(1)
    s = Server(inputQueue=inputQueue, outputQueue=outputQueue, host="localhost")
    print s
    s.run()

if __name__ == '__main__':
    
    # print possibleSigs((0,0,1,1))
    # 
    # sourcePath = '/Users/niko/arbeit/Data/SBMLmodels/iAF1260/biologValidatedSourcesWithElementalComposition.tsv'
    # sourceTable = readSourcesTable(sourcePath)
    # carbonSources = [elem[0] for elem in sourceTable if elem[1] == 1]
    # print carbonSources
    # sources2sig = dict([(row[0], tuple(row[1:]))for row in sourceTable])
    
    # testCombinatoricalSources(sourceTable, carbonSources)
    
    # combSources = list(combinatoricalSources(sourceTable, carbonSources))
    # print list(generateCombinatoricalMedia(20., combSources))[0:10]
    
    # include = dict([['R("R_ATPM")', [8.39, 8.39]],
    #                 ['R("Mo2b_Transp")', [0, 18.5]],
    #                 ['R("Mco2b_Transp")', [-1000, 1000]],
    #                 ['R("Mh2ob_Transp")', [-1000, 1000]],
    #                 ['R("Mhb_Transp")', [-1000, 1000]],
    #                 ['R("Mna1b_Transp")', [-1000, 1000]],
    #                 ['R("Mkb_Transp")', [-1000, 1000]],
    #                 ['R("Mca2b_Transp")', [-1000, 1000]],
    #                 ['R("Mcu2b_Transp")', [-1000, 1000]],
    #                 ['R("Mmg2b_Transp")', [-1000, 1000]],
    #                 ['R("Mzn2b_Transp")', [-1000, 1000]],
    #                 ['R("Mmobdb_Transp")', [-1000, 1000]],
    #                 ['R("Mfe2b_Transp")', [-1000, 1000]],
    #                 ['R("Mfe3b_Transp")', [-1000, 1000]],
    #                 ['R("Mcobalt2b_Transp")', [-1000, 1000]],
    #                 ['R("Mmn2b_Transp")', [-1000, 1000]],
    #                 ['R("Mclb_Transp")', [-1000, 1000]],
    #                 ['R("R_CAT")', [0, 0]],
    #                 ['R("R_SPODM")', [0, 0]],
    #                 ['R("R_SPODMpp")', [0, 0]],
    #                 ['R("R_FHL")', [0, 0]]])
    
    # gen = generateCombinatoricalMedia(20., combSources)
    # for i in range(10):
    #     medium = gen.next()
    #     print medium
    #     solveMedium('../models/iAF1260templateMinMax.lp', medium=medium, include=include, objective='R("R_Ec_biomass_iAF1260_core_59p81M")')


    try:
        sys.argv[1]
    except IndexError:
        sys.argv.append('server')
        sys.argv.append('test.h5')

    usage = """Usage:
python RandomMediaSimulations.py standalone storagefile configfile --> standalone mode
python RandomMediaSimulations.py server storagefile configfile --> server mode
python RandomMediaSimulations.py client serverip --> client mode"""
    try:
        if sys.argv[1] == 'standalone':
            basicFunctionality(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == 'server':
            server(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == 'client':
            client(sys.argv[2])
        else:
            print usage
    except IndexError:
        print usage
