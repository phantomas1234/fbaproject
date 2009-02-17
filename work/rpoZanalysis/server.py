#!/usr/bin/env python
# encoding: utf-8
"""
server.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import ifba.distFBA.networking as nw
import random
import Queue
import threading
import time
import pickle
import pprint
from ifba.combinatorics.combinatorics import SetCombine

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
    
def dataDelegateStub(text):
    scrambled = ''.join(random.sample(text, len(text)))
    return scrambled
    
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
    def __init__(self, queue):
        super(OutputClient, self).__init__()
        self.queue = queue
    
    def run(self):
        while True:
            print "Getting Data from Queue", self.queue
            print "The OutputQueue is that big: ", self.queue.qsize()
            data = self.queue.get()
            # pprint.pprint(data)
            file = open('log.txt', 'a')
            file.write(','.join([str(i) for i in data[1]]) + '\n')
            # pickle.dump(data, file)
            # for elem in data[1]:
            #     print elem
            #     # tmp = data[1][elem].tolist()
            #     print len(tmp)
            #     tmp2 = [str(i) for i in tmp]
            #     print tmp2
            #     file.write(', '.join(tmp2) + '\n')
            file.close

if __name__ == '__main__':
    generator = SetCombine(range(1,10),range(20,30)).generate()
    for i in generator:
        print i
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    t1 = InputClient(inputQueue, generator.next)
    t1.start()
    t2 = OutputClient(outputQueue)
    t2.start()
    s = nw.Server(inputQueue=inputQueue, outputQueue=outputQueue)
    print s
    s.run()    