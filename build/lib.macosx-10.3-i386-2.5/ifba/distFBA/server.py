#!/usr/bin/env python
# encoding: utf-8
"""
server.py

Created by Nikolaus Sonnenschein on 2008-01-22.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import ifba.distFBA.networking as nw
import random
import Queue
import threading
import time

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
            print "The queue is that big: ", self.queue.qsize()
            data = self.queue.put(self.func())
        
class OutputClient(threading.Thread):
    def __init__(self, queue):
        super(OutputClient, self).__init__()
        self.queue = queue
    
    def run(self):
        while True:
            print "Getting Data from Queue", self.queue
            print "The queue is that big: ", self.queue.qsize()
            data = self.queue.get()
            print data


if __name__ == '__main__':
    dataDelegateStub(delegateStub2())
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    t1 = InputClient(inputQueue, delegateStub2)
    t1.start()
    time.sleep(1)
    print "sadflkjas;dfkjsa;dlkfj;sladkjf"
    t2 = OutputClient(outputQueue)
    t2.start()
    s = nw.Server(inputQueue=inputQueue, outputQueue=outputQueue, port=50000)
    print s
    s.run()    