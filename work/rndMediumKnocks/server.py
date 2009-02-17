#!/usr/bin/env python
# encoding: utf-8
"""
server.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import ifba.distFBA.networking as nw
import sys
import random
import Queue
import threading
import time
import pickle
import pprint
from ifba.general.util import sumDicts, filterDict, dict2mathematica

def delegateStub2():
    string = """Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do
eiusmod tempor incididunt ut labore et dolore magna
aliqua. Ut enim ad minim veniam, quis nostrud exercitation
ullamco laboris nisi ut aliquip ex ea commodo consequat.
Duis aute irure dolor in reprehenderit in voluptate velit esse
cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat
cupidatat non proident, sunt in culpa qui officia deserunt mollit
anim id est laborum."""
    return string.swapcase()

def subMitThreshold():
    return sys.argv[2]

def dataDelegateStub(text):
    scrambled = ''.join(random.sample(text, len(text)))
    return scrambled
    
def dict2tsv(d):
    s = ''
    for key, value in d.items():
        s += str(key).replace('(', '[').replace(')', ']') + "\t" + repr(value) + "\n"
    return s

class InputClient(threading.Thread):
    def __init__(self, queue, func):
        super(InputClient, self).__init__()
        self.queue = queue
        self.func = func
    
    def run(self):
        while True:
            print "Putting data into queue", self.queue
            print "The InputQueue is that big: ", self.queue.qsize()
            data = self.queue.put(self.func())

class OutputClient(threading.Thread):
    def __init__(self, queue, pathPrefix):
        super(OutputClient, self).__init__()
        self.queue = queue
        self.prefix = pathPrefix
        self.activityDict = dict()
        self.koDict = dict()
    
    def run(self):
        counter = 0
        while True:
            counter = counter + 1
            print "Getting Data from Queue", self.queue
            print "The OutputQueue is that big: ", self.queue.qsize()
            data = self.queue.get()
            self.activityDict = sumDicts(self.activityDict, data[0])
            self.koDict = sumDicts(self.koDict, data[1])
            print self.koDict
            print self.activityDict
            path = self.prefix + "dict" + str(counter) + ".txt"
            print "writing to " + path
            open(path, 'w').write(dict2tsv(self.activityDict) + "\n" + dict2tsv(self.koDict) + "\n" + dict2tsv(data[2]))


if __name__ == '__main__':
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    t1 = InputClient(inputQueue, delegateStub2)
    t1.start()
    time.sleep(1)
    t2 = OutputClient(outputQueue, sys.argv[1])
    t2.start()
    s = nw.Server(inputQueue=inputQueue, outputQueue=outputQueue, host="localhost")
    print s
    s.run()