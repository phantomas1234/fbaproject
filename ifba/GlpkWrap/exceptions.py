#!/usr/bin/env python
# encoding: utf-8
"""
excpecptions.py

Created by Nikolaus Sonnenschein on 2009-12-04.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

class SolutionUnbounded(Exception):
    """docstring for SolutionUnbounded"""
    pass

class SolutionInfeasible(Exception):
    """docstring for SolutionUnbounded"""
    pass

class NoFeasibleSolution(Exception):
    """docstring for SolutionUnbounded"""
    pass

class SolutionUndefined(Exception):
    """docstring for SolutionUnbounded"""
    pass
