#!/usr/bin/env python
# encoding: utf-8
"""
installMaria.py

Created by Nikolaus Sonnenschein on 2008-01-16.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os


def main():
    os.system('cd glpk-4.36; ./configure; make; make install; make clean; \
    make distclean')
    os.system('python setup.py install')
    os.system('cd build; rm -R *')

if __name__ == '__main__':
    main()

