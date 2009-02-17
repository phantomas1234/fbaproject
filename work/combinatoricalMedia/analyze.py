#!/usr/bin/env python
# encoding: utf-8
"""
analyze.py

Created by Nikolaus Sonnenschein on 2008-06-04.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import pickle


def main():
    f = open('medKOmatrix.pic', 'rb')
    obj = pickle.load(f)
    f.close()
    # print obj[0][0]
    f = open('test.tsv', 'a')
    print len(obj)
    matrix = obj[0]
    for ar in matrix:
        f.write('\t'.join([str(i) for i in ar.tolist()]) + "\n")
    f.close()


if __name__ == '__main__':
    main()

