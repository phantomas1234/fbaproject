#!/usr/bin/env python
# encoding: utf-8
"""
util.py

Created by Nikolaus Sonnenschein on 2008-03-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import random
import string

def sumDicts(*dicts):
    """Sums up dicts.
    
    sumDicts({'Hanswurs':777777},
            {'Peter': 0.3}, 
            {}, 
            {'Hanswurs':-777775, 'Peter':0.9}) 
    ==> {'Hanswurs': 2, 'Peter': 1.2}
    """
    tmpDict = dict()
    keyList = list()
    keys = [e.keys() for e in dicts]
    for sublist in keys:
        keyList.extend(sublist)
    for key in set(keyList):
        values = []
        for elem in dicts:
            try:
                values.append(elem[key])
            except KeyError:
                values.append(0)
        tmpDict[key] = sum(values)
    return tmpDict
    
# def combineDicts(*dicts):
#     tmpDict = dict()
#     keyList = list()
#     keys = [e.keys() for e in dicts]
#     for sublist in keys:
#         keyList.extend(sublist)
#     for key in set(keyList):
#         values = []
#         for elem in dicts:
#             try:
#                 values.append(elem[key])
#             except KeyError:
#                 values.append(0)
#         tmpDict[key] = sum(values)
#     return tmpDict
    
def subDict(someDict, someKeys, default=None):
    return dict([(k, someDict.get(k, default)) for k in someKeys])

def filterDict(dict, threshold=0.):
    collect = [] 
    for key in dict.keys():
        if dict[key] > threshold:
            collect.append(key)
    return subDict(dict, collect)
            
def dict2mathematica(dict):
    """docstring for fluxDict2mathematica"""
    out = string.replace(str(dict), ':', ' ->')
    out = string.replace(out, "'", "")
    out = string.replace(out, "(", "[")
    out = string.replace(out, ")", "]")
    return out

if __name__ == '__main__':
    testDict = dict()
    for i, e in enumerate([random.uniform(0., 10.) for i in xrange(100)]):
        testDict[i] = e
    print "the TestDict", testDict
    print testDict[1]
    print sumDicts({'Hanswurs':777777}, {'Peter': 0.3}, {}, 
                    {'Hanswurs':-777775, 'Peter':0.9})                    
    print subDict(testDict, (1,4,7,77))
    
    print filterDict(testDict, 4.)
    print dict2mathematica(testDict)
    

# CLASS NAME: DLLInterface
#
# Author: Larry Bates (lbates@syscononline.com)
#
# Written: 12/09/2002
#
# Released under: GNU GENERAL PUBLIC LICENSE
#
#
class progressbarClass(object): 
    def __init__(self, finalcount, progresschar=None):
        import sys
        self.finalcount=finalcount
        self.blockcount=0
        #
        # See if caller passed me a character to use on the
        # progress bar (like "*").  If not use the block
        # character that makes it look like a real progress
        # bar.
        #
        if not progresschar: self.block=chr(178)
        else:                self.block=progresschar
        #
        # Get pointer to sys.stdout so I can use the write/flush
        # methods to display the progress bar.
        #
        self.f=sys.stdout
        #
        # If the final count is zero, don't start the progress gauge
        #
        if not self.finalcount : return
        self.f.write('\n------------------ % Progress -------------------1\n')
        self.f.write('    1    2    3    4    5    6    7    8    9    0\n')
        self.f.write('----0----0----0----0----0----0----0----0----0----0\n')
        return

    def progress(self, count):
        #
        # Make sure I don't try to go off the end (e.g. >100%)
        #
        count=min(count, self.finalcount)
        #
        # If finalcount is zero, I'm done
        #
        if self.finalcount:
            percentcomplete=int(round(100*count/self.finalcount))
            if percentcomplete < 1: percentcomplete=1
        else:
            percentcomplete=100
            
        #print "percentcomplete=",percentcomplete
        blockcount=int(percentcomplete/2)
        #print "blockcount=",blockcount
        if blockcount > self.blockcount:
            for i in range(self.blockcount,blockcount):
                self.f.write(self.block)
                self.f.flush()
                
        if percentcomplete == 100: self.f.write("\n")
        self.blockcount=blockcount
        return
    
def dict2tsv(d):
    s = ''
    for key, value in d.items():
        s += str(key).replace('(', '[').replace(')', ']') + "\t" + repr(value) + "\n"
    return s

def randomString(length):
    """Generates a random string of LENGTH length."""
    chars = string.letters + string.digits
    s = ""
    for i in random.sample(chars, length):
        s += i
    return s

def writeSimulationMetaData(path='.', notes=''):
    return os.system("svn info")


if __name__ == "__main__":
    print writeSimulationMetaData()
    # from time import sleep
    # pb=progressbarClass(8,"*")
    # count=0
    # while count<9:
    #     count+=1
    #     pb.progress(count)
    #     sleep(.2)
    # 
    # pb=progressbarClass(100, "*")
    # pb.progress(20)
    # sleep(0.2)
    # pb.progress(47)
    # sleep(0.2)
    # pb.progress(90)
    # sleep(0.2)
    # pb.progress(100)
    # print "testing 1:"
    # pb=progressbarClass(1, "*")
    # pb.progress(1)
    
