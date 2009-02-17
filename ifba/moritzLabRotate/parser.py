#!/usr/bin/env python
# encoding: utf-8
"""
parser.py

Created by Moritz Beber on 2008-10-13.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import re

def isAString(obj):
    return isinstance(obj, basestring)

def parseTwoColTSV2Dict(filename):
    if not os.path.isfile(filename):
        print "Error: File \"", filename, "\" does not exist!"
        return dict({})
    file = open(filename,"r")
    tmp = file.readlines()
    file.close()
    data = dict()
    for line in tmp:
        ltmp = line.split()
        if (len(ltmp) == 2):
            data[ltmp[0]] = ltmp[1]
    return dict(data)

def parseKEGGMetCentric2Dict(filename):
    if not os.path.isfile(filename):
        print "Error: File \"", filename, "\" does not exist!"
        return dict({})
    file = open(filename,"r")
    tmp = file.readlines()
    file.close()
    links = dict()
    for line in tmp:
        ltmp = line.split()
        if (len(ltmp) == 2):
            if ltmp[0] in links:
                links[ltmp[0]].append(ltmp[1])
            else:
                links[ltmp[0]] = [ltmp[1]]
    return dict(links)

def parseKEGGBipartite2Dict(filename):
    if not os.path.isfile(filename):
        print "Error: File \"", filename, "\" does not exist!"
        return dict({})
    file = open(filename,"r")
    tmp = file.readlines()
    file.close()
    links = dict()
    for line in tmp:
        ltmp = line.split()
        if (len(ltmp) == 2):
            if ltmp[0] in links:
                links[ltmp[0]].append(ltmp[1])
            else:
                links[ltmp[0]] = [ltmp[1]]
    return dict(links)

def parseMetaboliteKEGGMapping2Dict(filename):
    if not os.path.isfile(filename):
        print "Error: File \"",filename,"\" does not exist!"
        return dict({})
    file = open(filename,"r")
    tmp = file.readlines()
    file.close()
    abbr = re.compile("(\\w+)")
    kegg_id = re.compile("\\t(C\\d+)\\t")
    metabolites = dict()
    for line in tmp:
        mid = kegg_id.search(line)
        mabbr = abbr.match(line)
        if mid:
            metabolites["%s%s%s" % ('M', mabbr.group(1), 'c')] = mid.group(1)
            #print mabbr.group(1), mid.group(1)
        else:
            #print mabbr.group(1)
            pass
    return dict(metabolites)

if __name__ == '__main__':
    #test = parseKEGGMetCentric2Dict(sys.argv[1])
    #test = parseKEGGBipartite2Dict(sys.argv[1])
    #test = parseMetaboliteKEGGMapping2Dict(sys.argv[1])
    #test = parseTwoColTSV2Dict(sys.argv[1])
    #for key in test:
    #    print key, ":", test[key]
    pass

