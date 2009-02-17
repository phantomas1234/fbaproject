#!/usr/bin/env python
# encoding: utf-8
"""
setup.py

Created by Nikolaus Sonnenschein on 2008-01-07.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""


from distutils.core import setup, Extension

setup(
    name="ifba", 
    version="2.0",
    description="Flux Balance Analysis Software",
    author="Nikolaus Sonnenschein",
    author_email="niko.sonnenschein@googlemail.com",
    url="",
    # install only the swig library
    # packages=["ifba","ifba.glpki"],
    packages=["ifba","ifba.glpki", "ifba.GlpkWrap", "ifba.distributedFBA", 
    "ifba.general"],

    ext_modules=[Extension("_glpki", [("ifba/glpki/glpki.i",)],
                 include_dirs=["/Users/niko/arbeit/Software/glpk-4.32/include",
                 "/people/home/nsonnensch/downloads/glpk-4.32/include"],
                 library_dirs=["/usr/local/lib",
                 "/people/home/nsonnensch/software/lib"],
                 libraries=['glpk'],
                 extra_linker_args=['-arch i386'])
                 ]
    )

# include_dirs=["/home/engineer/software/glpk-4.32/include",
# "/people/home/nsonnensch/downloads/glpk-4.24/include",
# "/Users/niko/arbeit/Software/glpk-4.32/include"],