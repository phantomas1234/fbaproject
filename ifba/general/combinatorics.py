#!/usr/bin/env python
# encoding: utf-8
"""
combinatorics.py

Created by Nikolaus Sonnenschein on 2008-01-11.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

from bisect import *
from types import *
import random

class SetCombine(object):
    def __init__(self, list1, list2):
        # sorted set no. 1 which gets prepended to ...
        self.list1 = self.checkArg(list1)
        # sorted set no. 2
        self.list2 = self.checkArg(list2)
        # the last elements of subelements of no.1 (needed for bisect)
        self.list3 = tuple([i[-1] for i in self.list1])
        # the first elements of the subelements of no.2 (needed for bisect)
        self.list4 = tuple([i[0] for i in self.list2])
    
    def checkArg(self, l):
        # if type(l[0]) == IntType:
        #     return [(i) for i in l]
        # elif type(l[0]) == TupleType:
        #     return l
        # else:
        #     raise Exception, "The first Argument %s is not a list of tuples of \
        #     the form [1,2,3,4,5 ...] or [[1,2], [2,3] ...]" % l
        return l

    def generate(self):
        for i in range(0, len(self.list3)):
            pos = bisect(self.list4, self.list3[i])
            for e in range(pos, len(self.list4)):
                yield self.list1[i] + self.list2[e]

class SetCombineInfinite(SetCombine):
    """docstring for SetCombineInfinite"""
    def __init__(self, list1, list2):
        SetCombine.__init__(self, list1, list2)
        
    def generate(self):
        for i in range(0, len(self.list3)):
            pos = bisect(self.list4, self.list3[i])
            for e in range(pos, len(self.list4)):
                yield self.list1[i] + self.list2[e]


# def main():
#     list = [[i] for i in range(1,201)]
#     gen = SetCombine(list, list).generate()
#     # for i in gen:
#     #     print i
#     coll = [i for i in gen]
#     print len(coll)

def main2():
    list = [(i,) for i in range(200,202)]
    obj = SetCombine([(1,), (2,), (3,), (4,)], list).generate()
    for i in obj:
        print i

if __name__ == '__main__':
    main2()
    # g = SetCombineInfinite(range(1,20),range(10,30)).generate()
    # print list(g)

