#!/usr/bin/env python
# encoding: utf-8
"""
iAF1260simulations.py

Created by Nikolaus Sonnenschein on 2009-09-01.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import time
import Queue
import yaml
from ifba.storage.hdf5storage import SimulationDB, h5Container
from ifba.GlpkWrap.randomMedia import RandomMediaSimulations
from ifba.distributedFBA.networking import Server, Client
from ifba.distributedFBA.concurrency import stubInputClient, h5OutputClient


def readYamlConfig(path):
    """docstring for readYamlConfig"""
    config = yaml.load(open(path))
    return config

def generateRandomMediaObject(path=None, include=None, objective=None, minimizeRest=True):
    """Generate the RandomMediaSimulations object"""
    return RandomMediaSimulations(path, objective, include, '')

def generateStorageObject(path, config):
    """docstring for generateStorageObject"""
    randomSimulationsObj = generateRandomMediaObject(**config)
    return SimulationDB(h5Container(path, randomSimulationsObj.almaas.lp))

def basicFunctionality(outputfile, config):
    """Perform random media simulations and store them in a hdf5 file"""
    randomSimulationsObj = generateRandomMediaObject(**config)
    simulationStorage = generateStorageObject(outputfile, config)
    simulResult = randomSimulationsObj.run()
    simulationStorage.writeSimulationResult(simulResult)
    simulationStorage.close()

def client(serverip, configfile):
    """docstring for client"""
    config = readYamlConfig(configfile)
    randomSimulationsObj = generateRandomMediaObject(**config)
    counter = 0
    while True:
        counter = counter + 1
        print counter
        client = Client(task=randomSimulationsObj, host=serverip)
        client.run()

def server(outputfile='test.h5', config='parameters.yaml'):
    """Server"""
    config = readYamlConfig(config)
    def message2client():
        return "HooHoo"
    simulationStorage = generateStorageObject(outputfile, config)
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    t1 = stubInputClient(inputQueue, message2client)
    t1.start()
    time.sleep(1)
    t2 = h5OutputClient(outputQueue, simulationStorage)
    t2.start()
    s = Server(inputQueue=inputQueue, outputQueue=outputQueue, host="localhost")
    print s
    s.run()

if __name__ == '__main__':
    # basicFunctionality()
    # print readYamlConfig('parameters.yaml')
    try:
        sys.argv[1]
    except IndexError:
        sys.argv.append('server')
        sys.argv.append('test.h5')
    
    usage = 'Usage:\npython RandomMediaSimulations.py server storagefile configfile --> server mode\npython iAF1260simulations.py client serverip configfile --> client mode'
    try:
        if sys.argv[1] == 'server':
            server(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == 'client':
            client(sys.argv[2], sys.argv[3])
        else:
            print usage
    except IndexError:
        print usage
