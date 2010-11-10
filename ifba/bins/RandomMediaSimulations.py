#!/usr/bin/env python
# encoding: utf-8
"""
RandomMediaSimulaitons.py

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
from ifba.distributedFBA.concurrency import configInputClient, h5OutputClient
from ifba.GlpkWrap.util import ImportCplex
from ifba.GlpkWrap.metabolism import Metabolism


def readYamlConfig(path):
    """docstring for readYamlConfig"""
    config = yaml.load(open(path))
    return config

def generateRandomMediaObject(path=None, include=None, objective=None, minimizerest=True, optimizationRoutine='fba', koQ=True, descr=''):
    """Generate the RandomMediaSimulations object"""
    return RandomMediaSimulations(path, objective, include, descr, minimizerest, koQ=koQ, optimizationRoutine=optimizationRoutine)

def generateStorageObject(path, lp):
    """docstring for generateStorageObject"""
    return SimulationDB(h5Container(path, lp))

def basicFunctionality(outputfile, configpath, runs):
    """Perform random media simulations and store them in a hdf5 file"""
    config = readYamlConfig(configpath)
    descr = yaml.dump(config)
    print descr
    config['descr'] = descr
    randomSimulationsObj = generateRandomMediaObject(**config)
    simulationStorage = generateStorageObject(outputfile, randomSimulationsObj.almaas.lp)
    for i in range(runs):
        if i % 10 == 0:
            print i
        randomSimulationsObj = generateRandomMediaObject(**config)
        simulResult = randomSimulationsObj.run()
        del randomSimulationsObj
        simulationStorage.writeSimulationResult(simulResult)
    simulationStorage.close()

def client(serverip):
    """docstring for client"""
    counter = 0
    client = Client(task=generateRandomMediaObject, host=serverip)
    while True:
        counter = counter + 1
        print counter
        client.run()

def server(outputfile='test.h5', config='parameters.yaml'):
    """Server"""
    config = readYamlConfig(config)
    descr = yaml.dump(config)
    print descr
    config['descr'] = descr
    randomSimulationsObj = generateRandomMediaObject(**config)
    lp = Metabolism(ImportCplex(config['path']))
    simulationStorage = generateStorageObject(outputfile, lp)
    inputQueue = Queue.Queue(20)
    outputQueue = Queue.Queue(20)
    t1 = configInputClient(inputQueue, config)
    t1.start()
    time.sleep(1)
    t2 = h5OutputClient(outputQueue, simulationStorage)
    t2.start()
    s = Server(inputQueue=inputQueue, outputQueue=outputQueue, host="localhost")
    print s
    s.run()

if __name__ == '__main__':
    # print readYamlConfig('parameters.yaml')
    try:
        sys.argv[1]
    except IndexError:
        sys.argv.append('server')
        sys.argv.append('test.h5')
    
    usage = """Usage:
python RandomMediaSimulations.py standalone storagefile configfile runs --> standalone mode
python RandomMediaSimulations.py server storagefile configfile --> server mode
python RandomMediaSimulations.py client serverip --> client mode"""
    try:
        if sys.argv[1] == 'standalone':
            basicFunctionality(sys.argv[2], sys.argv[3], int(sys.argv[4]))
        elif sys.argv[1] == 'server':
            server(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == 'client':
            client(sys.argv[2])
        else:
            print usage
    except IndexError:
        print usage
