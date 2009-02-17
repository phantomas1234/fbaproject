#!/usr/bin/env python
# encoding: utf-8
"""
rpoZanalysis.py

Created by Nikolaus Sonnenschein on 2008-05-08.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import copy
import pprint

from ifba.combinatorics import combinatorics as comb
from ifba.GlpkWrap import metabolism, util, fluxdist
from ifba.blockedReactions.blockedReactions import analyseBlockedReactions
from ifba.general.util import progressbarClass
import pickle

def loadSource(path):
    """Loads a list of transportReactions. Format:
    R("Macgamb_Transp")
    R("Madnb_Transp")
    R("MalaDb_Transp")..."""
    file = open(path, 'r')
    sources = [line.strip() for line in file]
    file.close()
    return sources

def loadMinimalMed(path):
    """Loads a list of transportReactions. Format:
    R["R_ATPM"] -> {8.39, 8.39}
    R["Mo2b_Transp"] -> {0, 18.5}..."""
    file = open(path, 'r')
    content = eval(file.read())
    file.close()
    return dict(content)

def init(path):
    """A template gets initialized."""
    struct = util.ImportCplex(path)
    blub = metabolism.Metabolism(struct)
    blub.toggleVerbosity()
    return blub
    
class AnalyseKnockOutSpecificity(object):
    """Gets initialized with a lp model template and a set of deadly and
    viable sources (for the mutant)"""
    def __init__(self, lp, deadlySources, viableSources, sets2check=()):
        super(AnalyseKnockOutSpecificity, self).__init__()
        self.lp = lp
        self.deadlySources = deadlySources
        self.viableSources = viableSources
        if not sets2check:
            self.sets2check = self._unionOfActiveWtReactions()
        else:
            self.sets2check = sets2check
        print self.sets2check
        # print len(self.sets2check)
        
    def _unionOfActiveWtReactions(self):
        """Returns a set of reactions that have been active under the different
        circumstances provided by sources."""
        active = set()
        sources = self.deadlySources # + self.viableSources
        for source in sources:
            reacTmp = self._determineActiveReactions(source)
            active = active.union(reacTmp)
        return [(i,) for i in active]
    
    def _determineActiveReactions(self, carbon):
        """Determines the reactions which are active in the wt und the specified
        carbon source."""
        self.lp.modifyColumnBounds({carbon : (0, 20)})
        fluxDist = self.lp.fba()
        self.lp.initialize()
        print fluxDist.getActiveReactions()
        return fluxDist.getActiveColumns()

    def essentialQ(self, reactionSet):
        """Returns True if reaction is essential."""
        tmpBounds = dict()
        for r in reactionSet:
            tmpBounds[r] = (0., 0.)
        # print tmpBounds
        self.lp.modifyColumnBounds(tmpBounds)
        self.lp.simplex()
        self.lp.undo()
        objVal = self.lp.getObjVal()
        # print objVal
        if objVal <= 0.005:
            return True
        else:
            return False

    def checkReaction(self, carbon, reaction, spec=""):
        """docstring for checkReaction"""
        # print counter
        # print spec + carbon, "gets used as the carbon source!"
        self.lp.modifyColumnBounds({carbon : (0, 20)})
        bool = self.essentialQ(reaction)
        self.lp.initialize()
        return bool

    def checkKOs(self):
        # pb=progressbarClass(len(self.sets2check),"*")
        essentials = set()
        unessentials = set()
        fresh = set()
        counter = 0
        for reaction in self.sets2check:
            counter += 1
            print ("#" * 80) + "\n"
            print self.lp.translateIndices(reaction), "gets knocked out2!" + "\n"
            flag = False
            for carbon in self.deadlySources:
                print carbon, "----> deadly"
                if self.checkReaction(carbon, reaction):
                    unessentials.discard(reaction)
                    print "Essential! Checking viable carbon sources!"
                    for carbon2 in self.viableSources:
                        print "\t", carbon2, "----> viable"
                        if not self.checkReaction(carbon2, reaction):
                            essentials.add(reaction)
                        else:
                            print "\tThe reaction was essential under viable conditions. It gets discarded"
                            essentials.discard(reaction)
                            flag = True
                            break
                    break
                else:
                    print "Not Essential! Next deadly Carbon source!"
                    unessentials.add(reaction)
                    flu = set([(i,) for i in fluxdist.FluxDist(self.lp).getActiveColumns()])
                    diff = flu.difference(set(self.sets2check).union(fresh))
                    print "Fresh reactions:\n", diff
                    fresh = fresh.union(diff)
                    pass
                if flag:
                    break
                # pb.progress(counter)
        return essentials, unessentials, fresh

    def checkKOs2(self):
        # pb=progressbarClass(len(self.sets2check),"*")
        essentials = set()
        unessentials = set()
        counter = 0
        for reaction in self.sets2check:
            counter += 1
            print ("#" * 80) + "\n"
            print self.lp.translateIndices(reaction), "gets knocked out2!" + "\n"
            flag = False
            for carbon in self.deadlySources:
                print carbon, "----> deadly"
                if self.checkReaction(carbon, reaction):
                    unessentials.discard(reaction)
                    print "Essential! Checking viable carbon sources!"
                    for carbon2 in self.viableSources:
                        print "\t", carbon2, "----> viable"
                        if not self.checkReaction(carbon2, reaction):
                            essentials.add(reaction)
                        else:
                            print "\tThe reaction was essential under viable conditions. It gets discarded"
                            essentials.discard(reaction)
                            flag = True
                            break
                    break
                else:
                    print "Not Essential! Next deadly Carbon source!"
                    unessentials.add(reaction)
                    pass
                if flag:
                    break
                # pb.progress(counter)
        return essentials, unessentials


def loadTheModel():
    # loading the targets
    deadly = loadSource('./targetCarbonSources.txt')
    # loading the viable sources
    viable = loadSource('./viableCarbonSources.txt')
    # loading minimal medium conditions
    minMed = loadMinimalMed('./minimalMedium.txt')
    lp = init('../models/iAF1260template.lp')
    lp.modifyColumnBounds(minMed)
    # lp.simplex()
    # print lp.getObjVal()
    return copy.copy(lp), deadly, viable

def main():
    """Mainloop with class construct."""
    lpready, deadly, viable = loadTheModel()
    # seq = [1112, 3330, 2933, 3335, 2242, 1172, 2932, 1149, 2934, 1147, 3287, 1145, 2705, 2243, 2245, 2240, 2930, 3332, 2244]
    # analysis = AnalyseKnockOutSpecificity(lpready, deadly, viable, [(i,) for i in seq])
    analysis = AnalyseKnockOutSpecificity(lpready, deadly, viable)    
    essentials, unessentials, fresh = analysis.checkKOs()
    tmpFile = open('tmpFile.pc', 'w')
    pickle.dump((essentials, unessentials, fresh), tmpFile)
    tmpFile.close()

    print analysis.lp.translateIndices([i[0] for i in essentials])
    print analysis.lp.translateIndices([i[0] for i in unessentials])
    print analysis.lp.translateIndices([i[0] for i in fresh])
    print len(unessentials) + len(fresh)
        
def main2():
    lpready, deadly, viable = loadTheModel()
    tmpFile = open('tmpFile.pc', 'r')
    (essentials, unessentials, fresh) = pickle.load(tmpFile)
    tmpFile.close()
    
    print fresh

    new = list(unessentials.union(fresh))
    new.sort()
    print new
    doubles = comb.SetCombine(new, new).generate()
    analysis = AnalyseKnockOutSpecificity(lpready, deadly, viable, doubles)    
    essentials2, unessentials2 = analysis.checkKOs2()

    tmpFile = open('tmpFile2.pc', 'w')
    pickle.dump((essentials2, unessentials2), tmpFile)
    tmpFile.close()

def main3():
    lpready, deadly, viable = loadTheModel()


    tmpFile2 = open('tmpFile.pc', 'r')
    (essentials, unessentials, fresh) = pickle.load(tmpFile2)

    tmpFile = open('tmpFile2.pc', 'r')
    (essentials2, unessentials2) = pickle.load(tmpFile)


    tmpFile.close()
    
    new = list(unessentials)
    new.sort()

    new2 = list(unessentials2)
    new2.sort()
    print len(new)
    print len(new2)
    doubles = comb.SetCombine(new, new2).generate()
    print len(list(doubles))
    doubles = comb.SetCombine(new, new2).generate()
    # analysis = AnalyseKnockOutSpecificity(lpready, deadly, viable, doubles)    
    # essentials2, unessentials2 = analysis.checkKOs2()

    tmpFile = open('tmpFile3.pc', 'w')
    pickle.dump((essentials2, unessentials2), tmpFile)
    tmpFile.close()

if __name__ == '__main__':
    # seq = [(i,) for i in range(1,200)]
    # g = comb.SetCombine(seq, seq).generate()
    # for i in g:
    #     print i
    main()
