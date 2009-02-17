#!/usr/bin/env python
# encoding: utf-8
"""
test_ConstraintModelling.py

Created by Nikolaus Sonnenschein on 2008-02-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import unittest
from ifba.GlpkWrap import metabolism, util
import copy

class test_ConstraintModelling(unittest.TestCase):
    def setUp(self):
        self.lp = util.ImportCplex('test_data/model.lp')
        self.glp = metabolism.Metabolism(self.lp)
        
    def testCopy(self):
        """Tests if the the magic __copy__ methods is used correctly by the
        copy module. It tests specifically if the orginal and copied lp reside
        on different memory locations."""
        lpcopy = copy.copy(self.glp)
        self.assertNotEqual(lpcopy, self.glp)
        self.assertNotEqual(lpcopy.lp, self.glp.lp)
        self.assert_(isinstance(lpcopy, metabolism.Metabolism))
        
    def testPGKmutationIndices(self):
        """Tests if mutations by indices works."""
        glp = self.glp
        ind = glp.translateColumnNames(['R("R_PGK")', 'R("R_PGK_Rev")'])
        dict = {}
        for i in ind:
            dict[i] = (0., 0.)
        glp.modifyColumnBounds(dict)
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.5836778206)
        glp.undo()
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)

    def testPGKmutationColRemoval(self):
        """Tests if a specific mutation is correct."""
        glp = self.glp
        glp.deleteReactions(['R("R_PGK_Rev")'])
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.5836778206)
        glp.undo()
        glp.simplex()
        ind = glp.translateColumnNames(['R("R_PGK_Rev")'])
        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
        
    def testATPremoval(self):
        """Tests if a metabolite removal works correct."""
        glp = self.glp
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
        glp.deleteMetabolites(['Matpc', 'Madpc'])
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.)
        glp.undo()
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)

    def testPGKmutation(self):
        """Tests if a specific mutation is correct."""
        glp = self.glp
        glp.modifyColumnBounds({'R("R_PGK")':(0, 0), 'R("R_PGK_Rev")':(0, 0)})
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.5836778206)
        glp.undo()
        glp.undo()
        glp.simplex()
        self.assertAlmostEqual(glp.getObjVal(), 0.9259122)
    
    def testWrongTypeException(self):
        """Tests if a wrong formated bounds dictionary raises an exception."""
        glp = self.glp
        dict = {0.2:(0, 0), 'R("R_PGK_Rev")':(0, 0)}
        self.assertRaises(TypeError, glp.modifyColumnBounds, dict)
        
    def testGetObjective(self):
        """docstring for testGetObjective"""
        obj = self.glp.getObjectiveFunction()
        self.assertEqual(obj, {'R("R_BiomassEcoli")' : 1.})

    def testSetObjective(self):
        """docstring for testGetObjective"""
        self.glp.setObjectiveFunction({'R("R_PGK")' : 1., 'R("R_PGK_Rev")' : 1.})
        reference = {'R("R_PGK")' : 1., 'R("R_PGK_Rev")' : 1.}
        self.assertEqual(self.glp.getObjectiveFunction(), reference)
        self.glp.initialize()
        self.assertEqual(self.glp.getObjectiveFunction(), {'R("R_BiomassEcoli")' : 1.})

    def testSetReactionObjective(self):
        """docstring for testGetObjective"""
        self.glp.setReactionObjective('R("R_PGK")')
        self.assertEqual(self.glp.getObjectiveFunction(), {'R("R_PGK")' : 1.})
        self.glp.initialize()
        self.assertEqual(self.glp.getObjectiveFunction(), {'R("R_BiomassEcoli")' : 1.})

    def testSetReactionObjectiveMinRest(self):
        """docstring for testGetObjective"""
        self.glp.setReactionObjectiveMinimizeRest('R("R_PGK")')
        reacsTmp = self.glp.getColumnIDs()
        reference = zip(reacsTmp, [-1.e-8 for i in range(0, len(reacsTmp))])
        referenceDict = dict(reference)
        referenceDict['R("R_PGK")'] = 1.
        self.assertEqual(self.glp.getObjectiveFunction(), referenceDict)
        self.glp.initialize()
        self.assertEqual(self.glp.getObjectiveFunction(), {'R("R_BiomassEcoli")' : 1.})

    
if __name__ == '__main__':
    # unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(test_ConstraintModelling)
    unittest.TextTestRunner(verbosity=6).run(suite)