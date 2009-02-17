#!/usr/bin/env python
# encoding: utf-8
"""
static_degree.py

Created by Moritz Beber on 2008-10-13.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import parser

def inDegree(dict_links, metabolite):
    degree = 0
    for key in dict_links:
        degree += dict_links[key].count(metabolite)
    return degree

def outDegree(dict_links, metabolite):
    if metabolite in dict_links:
        return len(dict_links[metabolite])
    else:
        return 0

def degreeEss(dict_links, dict_nodes, dict_ess, filename, elist, mlist):
    ilist = []
    for key in dict_ess:
        if key in dict_nodes:
            indeg = inDegree(dict_links, dict_nodes[key])
            outdeg = outDegree(dict_links, dict_nodes[key])
            if indeg + outdeg != 0:
                ss = "%s\t%s\t%s\t%s\n" % (key, str(indeg), str(outdeg), dict_ess[key])
                ilist.append(ss)
            else:
                if not key in elist:
                    elist.append(key)
        else:
            if not key in mlist:
                mlist.append(key)
    output = open(filename, "w")
    output.writelines(ilist)
    output.close()

def main(argv):
    dict_links = parser.parseKEGGMetCentric2Dict(argv[0])
    dict_nodes = parser.parseMetaboliteKEGGMapping2Dict(argv[1])
    files = os.listdir(argv[2])
    elist = []
    mlist = []
    for name in files:
        if name.endswith("metbs"):
            dict_ess = parser.parseTwoColTSV2Dict(os.path.join(argv[2],name))
            filename = os.path.join(argv[3], name.split('.')[0] + ".tsv")
            degreeEss(dict_links, dict_nodes, dict_ess, filename, elist, mlist)
    elist.sort()
    mlist.sort()
    message = "List of metabolites that have no KEGG identifier: %s\nList of metabolites that cannot be mapped onto the KEGG graph: %s" % (str(mlist),str(elist))
    return message

if __name__ == '__main__':
    if not len(sys.argv) == 5:
        usage = """Usage:
python static_degree.py [metabolite centric network file : string] [metabolite mapping file : string] [source directory for essentiality files : string] [target directory for output : string]"""
        print usage
    else:
        sys.exit(main(sys.argv[1:]))

