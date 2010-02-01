#!/usr/bin/env python
# encoding: utf-8
"""
ConstraintModelling.py

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
    
    # def __str__(self):
    #     """Print some information regarding the lp model."""
    #     super(Metabolism, self).__init__(self)
    
    def __copy__(self):
        return Metabolism(super(Metabolism, self).__copy__().lp)
    
    def getReactions(self):
        """Returns a list of all reaction ids in the model."""
        tmpSet = set(self.getColumnIDs())
        return list(tmpSet.difference(set(self.getTransporters())))
    
    def getMetabolites(self):
        """Returns a list of all metabolite ids in the model."""
        return self.getRowIDs()
        # num = self.getNumRows()
        # return list([glpki.glp_get_row_name(self.lp, i) for i in range(1, num+1)])
        
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
        """docstring for getColumnBounds"""
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
    
    def setReactionObjectiveMinimizeRest(self, reaction, coeff=1.):
        """docstring for setReactionAsObjective"""
        rest = self.getColumnIDs()
        reactionDict = dict()
        for reac in rest:
            reactionDict[reac] = -1.e-8
        reactionDict[reaction] = coeff 
        self.setObjectiveFunction(reactionDict)
    
    def fba(self, method='simplex'):
        """Solve the Flux Balance Model and return a flux distribution
        instance."""
        eval("self."+method+"()")
        return ifba.GlpkWrap.fluxdist.FluxDist(self)
    
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
    
    def getFluxDict(self):
        return dict([(r, self.primalValues()[i]) for i, r in enumerate(self.getColumnIDs())])
    
    def getShadowPriceDict(self):
        dualVal = self.dualValues()
        return dict([(r, dualVal[i]) for i, r in enumerate(self.getColumnIDs())])
        
    def addMetaboliteDrains(self, metabolites):
        """Frees the row boundaries for metabolites"""
        boundDict = dict()
        for met in metabolites:
            # boundDict[met] = ('-inf', 'inf')
            boundDict[met] = (0., 'inf')
        self.modifyRowBounds(boundDict)
    


if __name__ == '__main__':
    def main():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/iAF1260template.lp'))
        print lp.getObjVal()
        lp.simplex()
        print lp.getObjVal()
        # print lp.getReactions()
        # print lp.getTransporters()
        transp = lp.getTransporters()
        print lp.translateColumnNames(transp)
        lp.modifyColumnBounds({'R("R_PGK")':(0, 20), 'R("R_ENO")':(0, 20)})
        print lp.history
    
    def main2():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        lp.simplex()
        print lp.getObjVal()
        lp.modifyColumnBounds({'R("R_PGK")':(0, 0), 'R("R_PGK_Rev")':(0, 0)})
        print lp.history
        lp.simplex()
        print lp.getObjVal()
        lp.undo()
        lp.simplex()
        print lp.getObjVal()
    
    def main3():
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print lp.translateColumnNames(['R("R_PGK")', 'R("R_PGK_Rev")'])
    
    def main4():
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        lp2 = copy.copy(lp)
        print lp.fba()
    
    def main5():
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        lp.setReactionObjective('R("R_BiomassEcoli")')
        print lp.getObjectiveFunction()
        
    def main6():
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print lp.getColumnBounds()
    
    def main7():
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        lp.simplexOLD()
        print lp.getObjVal()
        print lp.getReactions()
        newBounds = dict()
        # for i in range
        # lp.modifyColumnBounds({'R("R_PGK")':(0, 0), 'R("R_PGK_Rev")':(0, 0)})
        # lp.simplexOLD()
        # print lp.getObjVal()
        # lp.simplex()
        # print lp.getObjVal()
        # lp.undo()
        # lp.simplex()
        # print lp.getObjVal()
        
    def main8():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print lp.getObjVal()
        f1 = lp.fba()
        f11 = [i[1] for i in f1.getActiveFluxDist()]
        print lp.getObjVal()
        print lp.getObjective()
        lp.setReactionObjectiveMinimizeRest('R("R_BiomassEcoli")')
        print lp.getObjective()
        f2 = lp.fba()
        f22 = [i[1] for i in f1.getActiveFluxDist()]
        print lp.getObjVal()
        print [i - e for (i,e) in zip(f11,f22)]
        util.WriteCplex(lp, 'debug.lp')
    
    def main9():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print lp.getNumRows()
        # print lp.getMetabolites()
        lp.deleteReactions(['R("R_PGK")', 'R("R_PGK_Rev")'])
        lp.simplex()
        print lp.getObjVal()
        lp.undo()
        lp.simplex()
        print lp.getObjVal()

    def testAddColumn():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print lp.getNumRows()
        lp.addColumns({'F':(0.,999999.,{2:1.,10:-1.})})
        print lp.getColumnCoef(lp.getNumCols())
        lp.undo()
        print lp.getColumnCoef(lp.getNumCols())
        lp.deleteReactions(['R("R_PGK")', 'R("R_PGK_Rev")'])
        lp.simplex()
        print lp.getObjVal()
        lp.undo()

    def testDeleteColumn():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print glp_get_col_name(lp.lp, 25)
        print lp.getColumnCoef(25)
        print glpk.sparseList(lp.getColumnCoef(25))
        print glpki.glp_get_col_name(lp.lp, 26)
        print glpk.sparseList(lp.getColumnCoef(26))
        lp.deleteColumns([25])
        print glpki.glp_get_col_name(lp.lp, 25)
        print glpk.sparseList(lp.getColumnCoef(25))
        lp.undo()
        print glpki.glp_get_col_name(lp.lp, 25)
        print glpk.sparseList(lp.getColumnCoef(25))
        print glpki.glp_get_col_name(lp.lp, 1473)
        print glpk.sparseList(lp.getColumnCoef(1473))

    def testDeleteRow():
        """docstring for main"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        print glpki.glp_get_row_name(lp.lp, 25)
        print lp.getRowCoef(25)
        print glpk.sparseList(lp.getRowCoef(25))
        print glpki.glp_get_row_name(lp.lp, 26)
        print glpk.sparseList(lp.getRowCoef(26))
        lp.deleteRows([25])
        print glpki.glp_get_row_name(lp.lp, 25)
        print glpk.sparseList(lp.getRowCoef(25))
        lp.undo()
        print glpki.glp_get_row_name(lp.lp, 25)
        print glpk.sparseList(lp.getRowCoef(25))
        print glpki.glp_get_row_name(lp.lp, 1473)
        print glpk.sparseList(lp.getRowCoef(1473))
    
    def testFBAfunction():
        """docstring for testFBAfunction"""
        lp = Metabolism(util.ImportCplex('test_data/model.lp'))
        f2 = lp.fba('interiorPoint').getActiveFluxDist()
        print lp
        lp.getObjVal()
        # print f2
        # f1 = lp.fba().getActiveFluxDist()
        # print f1
        # print len(f1)
        # print len(f2)
        # print f1 == f2

    # testFBAfunction()

    # testAddColumn()
    # testDeleteColumn()
    # testDeleteRow()
    
    lp = Metabolism(util.ImportCplex('test_data/model.lp'))
    print lp.translateRowNames(['Macg5sac'])[0]