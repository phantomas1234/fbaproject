#!/usr/bin/env python
# encoding: utf-8
"""
rpoZinit.py

Created by Nikolaus Sonnenschein on 2008-02-04.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import clamvParser

hosts = clamvParser.clamvFilter().filterNumProcesses(threshold=0)
hosts = [i[0] + '.clamv.iu-bremen.de' for i in hosts]
# try:
#     hosts.remove('tlabterm.clamv.iu-bremen.de')
hosts = hosts[-3::]

def main(hosts):
    # os.system('python ./rpoZserver.py > /dev/null &')
    print hosts
    for i in hosts:
        print i
        os.system('ssh nsonnensch@%s "source .bash_profile;cd /people/home/nsonnensch/ifba/distFBA/;python client.py magenta.iuhb02.iu-bremen.de ../GlpkWrap/test_data/model.lp > /dev/null &" > /dev/null &' % i)
    sys.exit(1)

if __name__ == '__main__':
    print 'used hosts:', hosts
    main(hosts)

