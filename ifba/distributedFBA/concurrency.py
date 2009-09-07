#!/usr/bin/env python
# encoding: utf-8
"""
concurrency.py

Created by Nikolaus Sonnenschein on 2009-09-04.
Copyright (c) 2009 . All rights reserved.
"""

import threading

class InputClient(threading.Thread):
    def __init__(self, queue):
        super(InputClient, self).__init__()
        self.queue = queue

    def run(self):
        raise NotImplemented, "You've instantiated OutputClient directly although it is an abstract base class. The run method has to implemented in a sub class."


class OutputClient(threading.Thread):
    def __init__(self, queue):
        super(OutputClient, self).__init__()
        self.queue = queue

    def __str__(self):
        print "Getting Data from Queue", self.queue
        print "The OutputQueue is that big: ", self.queue.qsize()
        return ""

    def run(self):
        raise NotImplemented, "You've instantiated OutputClient directly although it is an abstract base class. The run method has to implemented in a sub class."

class stubInputClient(InputClient):
    def __init__(self, queue, func):
        InputClient.__init__(self, queue)
        self.func = func

    def run(self):
        while True:
            print "Putting data into queue", self.queue
            print "The InputQueue is that big: ", self.queue.qsize()
            data = self.queue.put(self.func())


class h5OutputClient(OutputClient):
    def __init__(self, queue, simulationDB):
        OutputClient.__init__(self, queue)
        self.simulationDB = simulationDB

    def run(self):
        counter = 0
        while True:
            counter = counter + 1
            print "Getting Data from Queue", self.queue
            print "The OutputQueue is that big: ", self.queue.qsize()
            simulationResult = self.queue.get()
            self.simulationDB.writeSimulationResult(simulationResult)


if __name__ == '__main__':
	main()

