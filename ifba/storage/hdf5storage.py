#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Nikolaus Sonnenschein on 2009-07-03.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

from tables import *
from ifba.GlpkWrap.fluxdist import FBAsimulationResult

def h5Container(fileName, lp, title="", complevel=5, expectedrows=1000000):
    """Constructs a h5 file container that is suitable for the storage of FBAsimulationResult objects."""
    numReacs = lp.getNumCols()
    class FBAsimulations(IsDescription):
        fluxactivity = Float32Col(shape=(numReacs), pos=0)
        lowerBounds = Float32Col(shape=(numReacs), pos=1)
        upperBounds = Float32Col(shape=(numReacs), pos=2)
        objective = Float32Col(shape=(numReacs), pos=3)
        timeStamp = Float32Col(pos=4)
        modelPath = StringCol(100, pos=5)
    h5 = openFile(fileName, 'w')
    h5.createArray(where=h5.root, name='rxnMapping', object=lp.getColumnIDs(), title='Index - ReactionID mapping')
    filters = Filters(complevel=5, complib='zlib', fletcher32=True, shuffle=True)
    h5.createTable(where=h5.root, description=FBAsimulations, name='simulations', title='Table containing the activities and media conditions for FBA simulations.', filters=filters, expectedrows=expectedrows)
    h5.flush()
    print h5
    h5.close()
    return fileName
    
class SimulationDB(object):
    """docstring for SimulationDB"""
    def __init__(self, h5file):
        self.h5container = openFile(h5file, 'r+')
    
    def writeSimulationResult(self, simulationResult):
        "Takes a SimulationResult object and store it to the database."
        table = self.h5container.root.simulations
        row = table.row
        colNames = table.colnames
        for col in colNames:
            try:
                attr = getattr(simulationResult, col)
            except NameError:
                print "Missing of a mandatory column field in the provided simulation result"
                sys.exit(-1)
            row[col] = attr
        row.append()
        self.h5container.flush()
        
    def retrieveSimulationResultByID(self):
        """docstring for re"""
        # raise NotImplementedError, "Do it looser!"
        table = self.h5container.root.simulations
        print table[1]
        
    def close(self):
        """Close the database."""
        self.h5container.flush()
        self.h5container.close()
        print "bye bye"
        

if __name__ == '__main__':
    pass