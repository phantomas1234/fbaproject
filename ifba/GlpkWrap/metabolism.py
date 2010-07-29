#!/usr/bin/env python
# encoding: utf-8
"""
metabolism.py

Created by Nikolaus Sonnenschein on 2008-02-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import re
import types
import copy
from ifba.GlpkWrap import glpk, util
from ifba.glpki import glpki
import ifba.GlpkWrap.fluxdist

class Metabolism(glpk.glpk):
    def __init__(self, lp):
        super(Metabolism, self).__init__(lp)
    
    def __copy__(self):
        return Metabolism(super(Metabolism, self).__copy__().lp)
    
    def getReactions(self):
        """Returns a list of all reaction ids in the model."""
        tmpSet = set(self.getColumnIDs())
        return list(tmpSet.difference(set(self.getTransporters())))
    
    def getMetabolites(self):
        """Returns a list of all metabolite ids in the model."""
        return self.getRowIDs()
        
    def getSubstratesAndProducts(self, reaction):
        """Return a tuple containing the substrates and products of reaction"""
        colIndex = self.translateColumnNames([reaction])[0]
        coeff = self.getColumnCoef(colIndex)
        substratesIndices = list()
        productsIndices = list()
        for i, elem in coeff.items():
            if elem > 0.:
                productsIndices.append(i)
            elif elem < 0.:
                substratesIndices.append(i)
        substratesIndices.sort()
        productsIndices.sort()
        return (tuple(self.translateRowIndices(substratesIndices)), tuple(self.translateRowIndices(productsIndices)))
    
    def addTransporters(self, rowIDs, postfix="_Transp"):
        newColumns = dict()
        for r in rowIDs:
            newColumns[r+postfix] = ('-inf', 'inf', {self.translateRowNames([r])[0]:1.})
        self.addColumns(newColumns)
        
    def addArtificalReaction(self, rowIDs, name="Biomass"):
        newColumns = dict()
        coeff = dict()
        for r in rowIDs:
            coeff[self.translateRowNames([r])[0]] = -1.
        newColumns[name] = (0., 1000, coeff) 
        self.addColumns(newColumns)
    
    def getTransporters(self, postfix='_Transp'):
        """Returns a list of all available transporters in the model. The
        postfix can be changed if transporters are stigmatized by a different
        signature."""
        def revMatch(str):
            patt = re.compile(".+" + postfix)
            if patt.search(str):
                return False
            else:
                return True
        return list([elem for elem in self.getColumnIDs() if not revMatch(elem)])
    
    def modifyColumnBounds(self, boundDict):
        """Adds dictionary of bounds, e.g. {R("R_PGK"):(0, 20),
        R("R_ENO"):(0, 20)}, to lp. Uses the memory of lp."""
        convDict = {}
        for i in boundDict:
            if type(i) == types.IntType:
                convDict[i] = boundDict[i]
            elif type(i) == types.StringType:
                convDict[self.translateColumnNames([i])[0]] = boundDict[i]
            else:
                raise TypeError, "Only integers or string identifiers are allowed"
        super(Metabolism, self).modifyColumnBounds(convDict)
    
    def modifyRowBounds(self, boundDict):
        """Adds dictionary of bounds, e.g. {'Matpc':(0, 20),
        'Madpc':(0, 20)}, to lp. Uses the memory of lp."""
        convDict = {}
        for i in boundDict:
            if type(i) == types.IntType:
                convDict[i] = boundDict[i]
            elif type(i) == types.StringType:
                convDict[self.translateRowNames([i])[0]] = boundDict[i]
            else:
                raise TypeError, "Only integers or string identifiers are allowed"
        super(Metabolism, self).modifyRowBounds(convDict)
    
    def getColumnBounds(self):
        """docstring for getColumnBounds"""
        bDict = super(Metabolism, self).getColumnBounds()
        return dict(zip(self.translateColumnIndices(bDict.keys()), bDict.values()))

    def getRowBounds(self):
        """docstring for getRowBounds"""
        bDict = super(Metabolism, self).getRowBounds()
        return dict(zip(self.translateRowIndices(bDict.keys()), bDict.values()))
    
    def getObjectiveFunction(self):
        """Returns the current Objective function."""
        objList = super(Metabolism, self).getObjective()
        reacs = self.getColumnIDs()
        tmp = [(reacs[i], coef) for i, coef in enumerate(objList) if coef != 0.]
        return dict(tmp)
    
    def setObjectiveFunction(self, reactionDict):
        """Set new Objective."""
        indexDict = dict()
        for reac in reactionDict:
            index = self.translateColumnNames([reac])[0]
            indexDict[index] = reactionDict[reac]
        super(Metabolism, self).setObjective(indexDict)
    
    def setReactionObjective(self, reaction, coeff=1.):
        """docstring for setReactionAsObjective"""
        self.setObjectiveFunction({reaction : coeff})
    
    def setReactionObjectiveMinimizeRest(self, reaction, coeff=1., factor=-1.e-8):
        """
        Sets a specified reaction as the new objective (maximization) and
        sets for all other reactions a negative factor coefficient, thus
        minimizing their usage. Transporters, which can have negative fluxes
        are not minimized.
        """
        rest = self.getColumnIDs()
        t = self.getTransporters()
        reactionDict = dict()
        for reac in rest:
            if reac in t:
                reactionDict[reac] = 0
            else:
                reactionDict[reac] = factor
        reactionDict[reaction] = coeff 
        self.setObjectiveFunction(reactionDict)

    def fba(self, method='simplex'):
        """Solve the Flux Balance Model and return a flux distribution
        instance."""
        eval("self."+method+"()")
        return ifba.GlpkWrap.fluxdist.FluxDist(self)
    
    def pFBA(self, obj, reactions=None, factor=-1.):
        """
        As implemented in:
        Lewis et al. Omic data from evolved E. coli are consistent with 
        computed optimal growth from genome-scale models. Mol Syst Biol (2010) 
        vol. 6 (1) pp.
        """
        self.setReactionObjective(obj, coeff=1.)
        mx = self.fba()[obj]
        self.modifyColumnBounds({obj:(mx, "inf")})
        objDict = dict()
        if not reactions == None:
            t = self.getTransporters()
            for rxn in reactions:
                if reac not in t:
                    objDict[rxn] = -1.
        else:
            for rxn in self.getReactions():
                objDict[rxn] = -1.
        self.setObjectiveFunction(objDict)
        # print self.getObjectiveFunction()
        return self.fba()
        
    def fluxSumAnalysis(self):
        """
        As implemented in:
        Chung und Lee. Flux-sum analysis: a metabolite-centric approach for 
        understanding the metabolic network. BMC Syst Biol (2009) vol. 3 pp. 117
        """
        pass

    def deleteReactions(self, listOfReactions):
        "Takes a list of reactions and constrains them to no flux."
        reactions = dict()
        for reaction in listOfReactions:
            reactions[reaction] = (0,0)
        self.modifyColumnBounds(reactions)
    
    def deleteReactionsFromStoich(self, listOfReactions):
        "Takes a list of reactions and removes them from the stoichiometry matrix."
        colIndexes = self.translateColumnNames(listOfReactions)
        self.deleteColumns(colIndexes)
    
    def deleteMetabolites(self, listOfMetabolites):
        "Takes a list of metabolites and removes them from the model."
        reactions = dict()
        for metabolite in listOfMetabolites:
            rowCoef = self.getRowCoef(self.translateRowNames([metabolite]).pop())
            for index, reaction in enumerate(rowCoef):
                if (reaction != 0.0):
                    reactions[self.translateColumnIndices([index + 1]).pop()] = (0,0)
        self.modifyColumnBounds(reactions)

    def deleteMetabolitesFromStoich(self, listOfMetabolites):
        "Takes a list of metabolites and removes them from the stoichiometry matrix."
        rowIndexes = self.translateRowNames(listOfMetabolites)
        self.deleteRows(rowIndexes)
    
    def getFluxDict(self):
        return dict([(r, self.primalValues()[i]) for i, r in enumerate(self.getColumnIDs())])
    
    def getShadowPriceDict(self):
        dualVal = self.dualValues()
        return dict([(r, dualVal[i]) for i, r in enumerate(self.getColumnIDs())])
        
    def addMetaboliteDrains(self, metabolites):
        """Frees the row boundaries for metabolites"""
        boundDict = dict()
        for met in metabolites:
            boundDict[met] = (0., 'inf')
        self.modifyRowBounds(boundDict)
    
    def freeMetabolites(self, metabolites):
        """Frees the row boundaries for metabolites"""
        boundDict = dict()
        for met in metabolites:
            boundDict[met] = ('-inf', 'inf')
        self.modifyRowBounds(boundDict)


if __name__ == '__main__':

    def main():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        obj = lp.getObjectiveFunction().keys()[0]
        print obj
        normalFluxDist = lp.fba()
        pFBAfluxDist = lp.pFBA(obj)
        stuff =  normalFluxDist.getActiveFluxDist()
        # print [i for i in stuff if not re.search('.*_Transp.*', i[0]) and not re.search('.*R_EX.*', i[0])]
        
        normalActiveReactions = normalFluxDist.getActiveReactions()
        print sum([j for i, j in normalFluxDist.getActiveFluxDist()])
        pFBAactiveReactions = pFBAfluxDist.getActiveReactions()
        print sum([j for i, j in pFBAfluxDist.getActiveFluxDist()])
        # print normalActiveReactions
        # print pFBAactiveReactions
        print set(normalActiveReactions) - set(pFBAactiveReactions)
        print set(pFBAactiveReactions) - set(normalActiveReactions)

    main()