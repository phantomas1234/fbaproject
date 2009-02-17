#!/usr/bin/env python
# encoding: utf-8
"""
build.py

Created by Nikolaus Sonnenschein on 2008-01-16.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os


def main():
    os.system('python setup.py build')
    # os.system('cd build; rm -R *')

if __name__ == '__main__':
    main()

