#!/usr/bin/env python
# encoding: utf-8
"""
knockOut_client.py

Created by Nikolaus Sonnenschein on 2008-01-22.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import ifba.distFBA.networking as nw

class delegateStub(object):
    def run(self):
        return """Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do 
eiusmod tempor incididunt ut labore et dolore magna
aliqua. Ut enim ad minim veniam, quis nostrud exercitation
ullamco laboris nisi ut aliquip ex ea commodo consequat.
Duis aute irure dolor in reprehenderit in voluptate velit esse
cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat
cupidatat non proident, sunt in culpa qui officia deserunt mollit
anim id est laborum."""


def main():
    print "Hello I am the client. I print small caps."
    print delegateStub()
    client = nw.Client(task=delegateStub(), host='212.201.48.108')
    print "And now my best friend the server."
    client.run()
    
def mainLoop():
    counter = 1
    while counter < 1000:
        main()
        print counter
        counter += 1
    

if __name__ == '__main__':
    mainLoop()