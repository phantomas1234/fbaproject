#!/usr/bin/env python
# encoding: utf-8
"""
glpk_Test.py

Created by Nikolaus Sonnenschein on 2008-02-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import unittest
import util
import glpk
import sys
import copy
import random
import pickle

class test_glpk(unittest.TestCase):
    def setUp(self):
        self.lp = util.ImportCplex('test_data/model.lp')
        self.glp = glpk.glpk(self.lp)
    
    def testSimplex(self):
        """Tests if the real glpk simplex function spits out the correct
        objective value for the iJR904 model under glucose minimal medium
        condition."""
        self.glp.simplex()
        obj = glpk.glp_get_obj_val(self.glp.lp)
        self.assertAlmostEqual(obj, 0.9259122)
    
    def testSimplex2(self):
        """test if the python glpk method works correct."""
        self.glp.simplex()
        obj = glpk.glp_get_obj_val(self.glp.lp)
        obj2 = self.glp.getObjVal()
        self.assertAlmostEqual(obj, obj2)
    
    def testModifyBounds(self):
        """Tests if the modifiyBounds method operates correctly and if the
        modifications are reversible using the glpk history."""
        self.glp.modifyColumnBounds({1 : (0, 0.5)})
        self.glp.simplex()
        obj = self.glp.getObjVal()
        self.assertAlmostEqual(obj, 0.5)
        self.glp.undo()
        self.glp.simplex()
        obj = self.glp.getObjVal()
        self.assertAlmostEqual(obj, 0.9259122)
        
    def testPickle(self):
        """Tests pickleability of glpk objects."""
        pickledLP = pickle.dumps(self.glp)
        unpickledLP = pickle.loads(pickledLP)
        self.glp.simplex()
        unpickledLP.simplex()
        self.assertAlmostEqual(self.glp.getObjVal(), unpickledLP.getObjVal())
    
    def testCopy(self):
        """Tests if the the magic __copy__ methods is used correctly by the
        copy module. It tests specifically if the orginal and copied lp reside
        on different memory locations."""
        lpcopy = copy.copy(self.glp)
        self.assertNotEqual(lpcopy, self.glp)
        self.assertNotEqual(lpcopy.lp, self.glp.lp)
    
    def testPresolOnOff(self):
        """Tests if one can savely shutdown the lp presolver. There are
        floating point precision differences between GLP_ON and GLP_OFF. So
        assertAlmostEqual is used to check if they really are only
        max precision errors"""
        # num = self.glp.getNumCols() # TODO someday you'll have to change this
        num = 10
        def func(float):
            if float < 0.:
                return 0.
            else:
                return float
        for i in range(1, num + 1):
            self.glp.smcp.presolve = glpk.GLP_ON
            self.glp.modifyColumnBounds({i : (0., 0.)})
            self.glp.simplex()
            res_on = func(self.glp.getObjVal())
            self.glp.undo()
            
            self.glp.smcp.presolve = glpk.GLP_OFF
            self.glp.modifyColumnBounds({i : (0., 0.)})
            self.glp.simplex()
            res_off = func(self.glp.getObjVal())
            self.glp.undo()
            
            self.assertAlmostEqual(res_on, res_off)
    
    def testInitialize(self):
        """Tests if initialize returns the same clean lp object as a series
        of consecutive undo calls."""
        sample = random.sample(range(1, self.glp.getNumCols() + 1), 69)
        rndReal = random.uniform
        for i in sample:
            real1 = rndReal(20., 100.)
            self.glp.modifyColumnBounds({i:(0., real1)})
        self.glp.initialize()
        # self.glp.undo()
        # self.glp.undo()
        self.glp.simplex()
        obj = self.glp.getObjVal()
        self.assertAlmostEqual(obj, 0.9259122)
    
    def testObjectiveFunctionality(self):
        """Tests if the getter and setter function for the lp objective work
        and if they really are reversible."""
        initObjList = self.glp.getObjective()
        self.glp.setObjective({2 : 1.})
        self.assertNotEqual(self.glp.getObjective(), initObjList)
        self.glp.initialize()
        self.assertEqual(self.glp.getObjective(), initObjList)
    
    def testSetOptFlag(self):
        """Tests if the getter and setter functions for optimization direction
        work and if they are reversible."""
        flag = self.glp.getOptFlag()
        self.assertEqual(flag, 2)
        self.assertEqual(flag, glpk.GLP_MAX)
        self.glp.setOptFlag('MiN')
        flag = self.glp.getOptFlag()
        self.assertEqual(flag, glpk.GLP_MIN)
        self.glp.initialize()
        flag = self.glp.getOptFlag()
        self.assertEqual(flag, glpk.GLP_MAX)
    
    def testGetBounds(self):
        """Tests if the getColumnBounds method return a correct dictionary of column
        bounds."""
        # boundsDict = self.glp.getColumnBounds()
        # self.assertEqual(boundsDict)
        pass
    
    def testDeleteColumns(self):
        """Tests if a the specified column is deleted from the constraint
        matrix"""
        self.assertEqual(self.glp.getNumCols(), 1473)
        colCoef = self.glp.getColumnCoef(66)
        self.glp.deleteColumns([66])
        self.assertEqual(self.glp.getNumCols(), 1472)
        # now we check if this can be undone
        self.glp.undo()
        self.assertEqual(self.glp.getColumnCoef(1473), colCoef)
        self.assertEqual(self.glp.getNumCols(), 1473)
    
    def testDeleteRows(self):
        """Tests if a the specified row is deleted from the constraint
        matrix"""
        self.assertEqual(self.glp.getNumRows(), 904)
        rowCoef = self.glp.getRowCoef(800)
        self.glp.deleteRows([800])
        self.assertEqual(self.glp.getNumRows(), 903)
        # now we check if this can be undone
        self.glp.undo()
        self.assertEqual(self.glp.getNumRows(), 904)
        self.assertEqual(self.glp.getRowCoef(904), rowCoef)
    
    def testAddColumns(self):
        """Tests if a the specified column is appended to the constraint
        matrix"""
        self.assertEqual(self.glp.getNumCols(), 1473)
        newColumArray = self.glp.getColumnCoef(1)
        self.glp.addColumns({'R("R_HansWurs")': (0., 99999., newColumArray)})
        self.assertEqual(self.glp.getNumCols(), 1474)
        # now we check if this can be undone
        self.glp.undo()
        self.assertEqual(self.glp.getNumCols(), 1473)
        self.assertEqual(len(self.glp.history), 0)
    
    def testAddRows(self):
        """Tests if a the specified column is appended to the constraint
        matrix"""
        self.assertEqual(self.glp.getNumRows(), 904)
        newColumArray = self.glp.getRowCoef(1)
        self.glp.addRows({'Mwurstb': (0., 99999., newColumArray)})
        self.assertEqual(self.glp.getNumRows(), 905)
        # now we check if this can be undone
        self.glp.undo()
        self.assertEqual(self.glp.getNumRows(), 904)
        self.assertEqual(len(self.glp.history), 0)
    
    def testCheckIndexValidities(self):
        """Check if an IndexError is raised if an column or row index is out
        of range."""
        self.assertRaises(IndexError, self.glp._setColumnBound, 1474, 0., 0.)
        self.assertRaises(IndexError, self.glp._setColumnBound, 0, 0., 0.)        
        self.assertRaises(IndexError, self.glp._setRowBound, 905, 0., 0.)
        self.assertRaises(IndexError, self.glp._setRowBound, 0, 0., 0.)
        self.assertRaises(IndexError, self.glp.translateColumnIndices, [0])
        self.assertRaises(IndexError, self.glp.translateColumnIndices, [1474])
        self.assertRaises(IndexError, self.glp.translateRowIndices, [0])
        self.assertRaises(IndexError, self.glp.translateRowIndices, [905])

    def testTranslateRowColumnNames(self):
        """docstring for testTranslateRowColumnNames"""
        self.assertEqual(self.glp.translateColumnNames(['R("R_BiomassEcoli")', 
        'R("R_XYLI1_Rev")']), [1,1473])
        self.assertRaises(Exception, self.glp.translateColumnNames, ['R("R_Stub")'])
        self.assertEqual(self.glp.translateRowNames(['Matpc', 'MglcDb']), [223, 430])
        self.assertRaises(Exception, self.glp.translateRowNames, ['Mstubc'])
        
        
class test_sparseList(unittest.TestCase):
    def setUp(self):
        self.spL = glpk.sparseList([1, -3, 0, 0, -6., 0, 3, 0, 0, 0, 0., 3.])
        
    def testIt(self):
        self.assertEqual(self.spL, {1:1,2:-3,5:-6.,7:3,12:3.})
        
if __name__ == '__main__':
    # unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(test_glpk)
    unittest.TextTestRunner(verbosity=6).run(suite)
    suite = unittest.TestLoader().loadTestsFromTestCase(test_sparseList)
    unittest.TextTestRunner(verbosity=6).run(suite)