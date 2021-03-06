#!/usr/bin/env python
# encoding: utf-8
"""
server.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
Modified by Moritz Beber on 2008-11-26.
"""


import ifba.distFBA.networking as nw
import sys
import os
import gzip
import re
import random
import Queue
import threading
import time
import pickle
import pprint
from ifba.general.util import sumDicts, filterDict, dict2mathematica


def parseFluxFile(task):
    file = gzip.open(task,"r")
    tmp = file.read()
    file.close() 
    # extract only part before empty line
    pattern = re.compile('^\\n', re.M)
    tmp = re.split(pattern, tmp)
    # split single lines and transform braces
    tmp = tmp[0].splitlines() # use splitlines() rather than split('\n') to remove empty last line
    tmp = [s.replace('{','(').replace('}',')').split('\t') for s in tmp]
    # transforms string into proper number pair
    tmp = [(i[0], eval(i[1]))for i in tmp]
    return dict(tmp)


class InputClient(threading.Thread):
    def __init__(self, queue, tasks, directory, func):
        super(InputClient, self).__init__()
        self.queue = queue
        self.tasks = tasks
        self.directory = directory
        self.func = func
    
    def run(self):
        while len(self.tasks) > 0:
            print "The InputQueue has size:", self.queue.qsize()
            print "Putting data into queue", self.queue
            filename = self.tasks.pop()
            task = os.path.join(self.directory, filename)
            medium = self.func(task)
            name = filename.split('.')
            if os.path.isfile(task):
                self.queue.put((medium, name[0]))
            else:
                error = "Error: File \"%s\" not found!" % (task)
                print sys.stderr, error


class OutputClient(threading.Thread):
    def __init__(self, queue, directory):
        super(OutputClient, self).__init__()
        self.queue = queue
        self.directory = directory
    
    def run(self):
        counter = 0
        while True:
            counter += 1
            print "The OutputQueue is that big: ", self.queue.qsize()
            print "Getting Data from Queue", self.queue
            (rdict, mdict, name) = self.queue.get()
            filename = name + ".rxns"
            filename = os.path.join(self.directory, filename)
            outfile = open(filename, "w")
            for key in rdict.keys():
                s = "%s\t%s\n" % (key, str(float(rdict[key])))
                outfile.write(s)
            outfile.close()
            filename = name + ".metbs"
            filename = os.path.join(self.directory, filename)
            outfile = open(filename, "w")
            for key in mdict.keys():
                s = "%s\t%s\n" % (key, str(float(mdict[key])))
                outfile.write(s)
            outfile.close()


def main(argv):
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    # make list with all flux data files in directory given by argv
    files = os.listdir(argv[0])
    # make input client with job queue, files to work on, their directory, and the function to use
    t1 = InputClient(inputQueue, files, argv[0], parseFluxFile)
    t1.start()
    time.sleep(1)
    # start output client with job queue, output directory, info joining function
    t2 = OutputClient(outputQueue, argv[1])
    t2.start()
    s = nw.Server(inputQueue=inputQueue, outputQueue=outputQueue, host="localhost", port=50001)
    print s
    s.run()


if __name__ == '__main__':
    if (len(sys.argv) == 3):
        main(sys.argv[1:])
        print "Terminated normally."
        sys.exit(0)
    else:
        print >> sys.stderr, "Usage:\npython server.py [flux source directory] [data target directory]"
        sys.exit(2)

