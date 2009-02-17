#!/usr/bin/env python
# encoding: utf-8
"""
analysis.py

Created by Moritz Beber & Nikolaus Sonnenschein on 2008-09-26.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism, randomMedia, knockouts, fluxdist
import gzip
import re


def parseFluxFile(filename):
    """
    Reads the medium conditions from 'filename' which is to be in a specialised
    file format. Please check for file existance before calling.
    
    """
    file = gzip.open(filename,"r")
    tmp = file.read()
    file.close() 
    # extract only part before empty line
    pattern = re.compile('^\\n', re.M)
    tmp = re.split(pattern, tmp)
    # split single lines and transform braces
    tmp = tmp[0].splitlines() # use splitlines() rather than split('\n') to remove empty last line
    tmp = [s.replace('{','(').replace('}',')').split('\t') for s in tmp]
    # transforms string into proper number pair
    tmp = [(i[0], eval(i[1]))for i in tmp]
    return dict(tmp)
    

def predefMediumFBA(medium, name, target, lp):
    """
    Given the specific dictionary 'medium' this function performs a FBA on the
    iAF1260 model under these conditions. Also performs FBA on knock-outs.
    """
    lp.modifyColumnBounds(medium)
    lp.glpSimplex()
    fd = fluxdist.FluxDist(lp)
    wild_type = lp.getObjVal()
    if (wild_type == 0):
        error = "No biomass production in wild type under medium conditions \"%s\". Skipping..." % (name)
        print >> sys.stderr, error
    else:
        # knock-out each of the active reactions
        rxns = fd.getActiveReactions()
        name = name.split('.')
        filename = '%s%s%s' % (target, name[0], ".rxns")
        llist = []
        for rxn in rxns:
            #print rxn, 'is knocked out'
            lp.deleteReactions([rxn])
            try:
                lp.glpSimplex()
                if lp.getObjVal() < 0.0:
                    s = "%s\t%s\n" % (rxn, str(1.0))
                else:
                    s = "%s\t%s\n" % (rxn, str((wild_type - lp.getObjVal()) / wild_type))
                llist.append(s)
                lp.undo()
            except Exception, msg:
                print >> sys.stderr, msg
                print >> sys.stderr, "Knock-out of reaction", rxn, "has no feasable solution! Assuming 0 growth."
                s = "%s\t%s\n" % (rxn, str(1.0))
                llist.append(s)
                lp.undo()
        # write essentiality information to file
        output = open(filename, "w")
        output.writelines(llist)
        output.close()
        # knock-out each metabolite mapping of data onto active part needed later
        metbs = fd.getActiveMetabolites()
        filename = '%s%s%s' % (target, name[0], ".metbs")
        llist = []
        for metb in metbs:
            if (metb.endswith('c')):
                #print metb, 'is knocked out.'
                lp.deleteMetabolites([metb])
                try:
                    lp.glpSimplex()
                    if lp.getObjVal() < 0.0:
                        s = "%s\t%s\n" % (metb, str(1.0))
                    else:
                        s = "%s\t%s\n" % (metb, str((wild_type - lp.getObjVal()) / wild_type))
                    llist.append(s)
                    lp.undo()
                except Exception, msg:
                    print >> sys.stderr, msg
                    print >> sys.stderr, "Knock-out of metabolite", metb, "has no feasable solution! Assuming 0 growth."
                    s = "%s\t%s\n" % (metb, str(1.0))
                    llist.append(s)
                    lp.undo()
        # write essentiality information to file
        output = open(filename, "w")
        output.writelines(llist)
        output.close()
    # reset modifications for next medium conditions
    lp.initialize()

# requires a single directory with only flux information files
def importMedium(directory, target):
    """
    Loops through all files in 'directory' reading the media conditions and
    initiating essentiality analyses on each. Please have only appropriate files
    in the directory containing the medium information in distinct format.
    """
    lp = metabolism.Metabolism(util.ImportCplex('../../models/iAF1260template2.lp', terminal="OFF"))
    files = os.listdir(directory)
    for name in files:
        filename = os.path.join(directory, name)
        info = "Using flux information from: \"%s\"." % (name)
        print info
        # if the file exists it extracts the medium and starts analysis on it
        if os.path.isfile(filename):
            medium = parseFluxFile(filename)
            predefMediumFBA(medium, name, target, lp)
        else:
            error = "Error: File \"%s\" not found!" % (filename)
            print sys.stderr, error

if __name__ == '__main__':
    if (len(sys.argv) == 3):
        importMedium(sys.argv[1], sys.argv[2])
        print "Terminated normally."
        sys.exit(0)
    else:
        print >> sys.stderr, "Usage:\npython analysis.py [flux source directory] [data target directory]"
        sys.exit(2)

