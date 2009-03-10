#!/usr/bin/env python
# encoding: utf-8
"""
arndtColab.py

Created by Nikolaus Sonnenschein on 2009-03-10.
Copyright (c) 2009 Jacobs University Bremen. All rights reserved.
"""

import sys
import os
from ifba.ContextFBA import ContextFBA



if __name__ == '__main__':
    ContextFBA.cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/", \
    "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/eberResults2.tsv")

