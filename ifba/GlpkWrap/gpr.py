#!/usr/bin/env python
# encoding: utf-8
"""
gpr.py

This module and its class emulate the GPR (Gene-Protein-Reaction) associations
as represented in the BIGG database (http://bigg.ucsd.edu/)

Schellenberger et al. BiGG: a Biochemical Genetic and Genomic knowledgebase of
large scale metabolic reconstructions. BMC Bioinformatics (2010) vol. 11 (1)
pp. 213

Created by Nikolaus Sonnenschein on 2010-07-29.
Copyright (c) 2010 . All rights reserved.
"""

import networkx as nx

# =====================
# = Class definitions =
# =====================

class Entity(object):
    """docstring for Entity"""
    def __init__(self, id=None):
        super(Entity, self).__init__()
        self.id = id

class Locus(Entity):
    """d"""
    def __init__(self, id=None):
        super(Locus, self).__init__(id=id)

class Gene(Entity):
    """docstring for gene"""
    def __init__(self, id=None):
        super(Gene, self).__init__(id=id)

class Protein(Entity):
    """docstring for Protein"""
    def __init__(self, id=None):
        super(Protein, self).__init__(id=id)
        self.id = id

class Reaction(Entity):
    """docstring for Protein"""
    _instances = dict()
    def __init__(self, id=None):
        if not id in _instances:
            super(Reaction, self).__init__(id=id)
            self.id = id
            _instances[id] = self
        else:
            self.__dict__


class GPRmapping(object):
    """docstring for GPRmapping"""
    def __init__(self, arg):
        super(GPRmapping, self).__init__()
        self.arg = arg

# =====================
# = Utility Functions =
# =====================

def parseGpr(str):
    """docstring for parseGpr"""
    pass
