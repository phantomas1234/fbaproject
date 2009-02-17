#!/usr/bin/env python
# encoding: utf-8
"""
startClients.py

Created by Nikolaus Sonnenschein on 2008-10-09.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import urllib
import re
import string
import pprint
import os
import random

class clamvParser:
    """docstring for clamvParser"""
    def __init__(self):
        self.url = 'http://www.clamv.iu-bremen.de/CLAMV/Wizard/Accounting/html/CLAMVCluster.html'
        self.page = urllib.urlopen(self.url).read()
        regEx = re.compile('<! build Table 1 in the html file>.+<! end of table 1>', re.S)
        self.lines = string.split(regEx.search(self.page).group(), "\n")[7:-2]
        reg = re.compile('<td>(.+?)</td>')
        self.table = [reg.findall(elem)for elem in self.lines]
    
    def __str__(self):
        print self.table

class clamvFilter(clamvParser):
	def _sortBy_Load(self, x, y):
		if x[3] < y[3]:
			return 1
		elif x[3] == y[3]:
			return 0
		else:
			return -1
	
	def filterClusterLoad(self, threshold=0.):
		tab = filter(lambda x: float(x[-3]) <= threshold, self.table)
		tab.sort(self._sortBy_Load)
		return tab

	def filterNumProcesses(self, threshold=0):
		return filter(lambda x: float(x[-2]) <= threshold, self.table)

def helperFunc(table, oldOnesThresh=40):
    l = list()
    print table
    for i in table:
        tmp = i[1:]
        if i[0] == 'tlabterm':
            break
        tmp.insert(0, i[0][4:])
        l.append(tmp)
    return [i for i in l if int(i[0]) > oldOnesThresh]

def getHosts(numProc):
    clamv = clamvFilter()
    # print clamv
    # pprint.pprint(clamv.filterClusterLoad(1.))
    output = helperFunc(clamv.filterNumProcesses(numProc))
    hosts = ['tlab'+i[0]+'.clamv.iu-bremen.de' for i in output]
    return hosts
    
def startClient(host):
    # os.system('ssh nsonnensch@' + host + ' "cd Sandbox; ./remoteTest &"')
    # os.system('ssh nsonnensch@' + host + ' bash -c "Sandbox/remoteTest &" &')
    c = ' "source .bash_profile ; cd ./Sandbox/ifba/ifba/rndMediumKnocks/ ; python client.py 212.201.48.108 > /dev/null &"' 
    os.system('ssh nsonnensch@' + host + c + ' &')    


if __name__ == '__main__':
    hostsZero = set(getHosts(0))
    print hostsZero
    print len(hostsZero)
    hostsOne = set(random.sample(set(getHosts(1)).difference(hostsZero),10))
    print hostsOne.union(hostsZero)
    print len(hostsOne.union(hostsZero))
    superHosts = list(hostsOne.union(hostsZero))
    for host in superHosts:
        print 'starting client on ', host
        startClient(host)
    print "bye, bye my dear"
    exit()
