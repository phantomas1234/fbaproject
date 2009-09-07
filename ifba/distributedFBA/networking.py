#!/usr/bin/env python
# encoding: utf-8
"""
networking.py

Created by Nikolaus Sonnenschein on 2008-04-21.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import socket
import pickle
import threading
import select
import Queue, heapq, time
from ifba.general.util import sumDicts, filterDict, dict2mathematica
from optparse import OptionParser

parser = OptionParser()

import socket,struct,sys,time

class PriorityQueue(Queue.Queue):
    """From the Python Cookbook sec. edition"""
    def _init(self, maxsize):
        self.maxsize = maxsize
        self.queue = []
        
    def _qsize(self):
        """Return the number of items that are currently enqueued"""
        return len(self.queue)
    
    def _empty(self):
        """Check wether the queue is empty"""
        return not self.queue
        
    def _full(self):
        """docstring for _full"""
        return self.maxsize > 0 and len(self.queue) >= self.maxsize
        
    def _put(self, item):
        """Put a new item in the queue"""
        heapq.heappush(self.queue, item)
        
    def _get(self):
        """Get an item from the queue"""
        return heapq.heappop(self.queue)
        
    def put(self, item, priority=0, block=True, timeout=None):
        """docstring for put"""
        decorated_item = priority, time.time(), item
        Queue.Queue.put(self, decorated_item, block, timeout)
        
    def get(self, block=True, timeout=None):
        """docstring for get"""
        priority, time_posted, item = Queue.Queue.get(self, block, timeout)
        return item

class Networking(object):
    """Networking is an abstract base class for everything that has certain
    commonalities e.g. a hostname and a port etc.."""
    def __init__(self, host='localhost', port=50000, size=1042, 
    endMarker='!@#$^&*', *args, **kwargs):
        super(Networking, self).__init__(*args, **kwargs)
        self.host = host
        self.port = port
        self.size = size
        self.endMarker = endMarker
    
    def recv_end(self):
        total_data=[];data=''
        while True:
            data = self.sock.recv(self.size)
            if self.endMarker in data:
                total_data.append(data[:data.find(self.endMarker)])
                break
            total_data.append(data)
            if len(total_data)>1:
                #check if end_of_data was split
                last_pair=total_data[-2]+total_data[-1]
                if self.endMarker in last_pair:
                    total_data[-2]=last_pair[:last_pair.find(self.endMarker)]
                    total_data.pop()
                    break
        return ''.join(total_data)
        
    def recv_pickle(self):
        data = ''
        inp = True
        while inp:
            inp = self.sock.recv(self.size)
            data += inp
        return pickle.loads(data)

    def send_pickle(self, data):
        self.sock.sendall(pickle.dumps(data))
        
    def send_end(self, data):
        self.sock.sendall(data + self.endMarker)


class Transaction(Networking):
    """Transaction objects that can be passed to a thread. 
    
    They implement a run method that does everything the thread should do and 
    they include the information, which is necessary to all the networking 
    communication."""
    def __init__(self, sock=None, inputQueue=None, outputQueue=None, *args, **kwargs):
        super(Transaction, self).__init__(*args, **kwargs)
        self.sock = sock
        self.inputQueue = inputQueue
        self.outputQueue = outputQueue
    
    def run(self):
        """docstring for run"""
        print "sending to client", self.sock
        # self.send_end(self.inputQueue.get())
        # self.sock.sendall(self.inputQueue.get())
        self.send_pickle(self.inputQueue.get())
        self.sock.shutdown(socket.SHUT_WR)
        print "recieving from client", self.sock
        # fromClient = self.recv_end()
        fromClient = self.recv_pickle()
        self.outputQueue.put(fromClient)
        self.sock.close()


class Server(Networking):
    def __init__(self, inputQueue=None, outputQueue=None, *args, **kwargs):
        super(Server, self).__init__(*args, **kwargs)
        try:
            self.host = socket.gethostbyname(socket.gethostname())
        except:
            print "hostname not available"
            print "set host to 'localhost'"
            self.host = 'localhost'
        self.backlog = 5
        self.threads = []
        self.sock = None
        self.inputQueue = inputQueue
        self.outputQueue = outputQueue

    
    def __str__(self):
        return """Some information about the server:
hostname -> %s
port -> %s
number of active threads -> %s
""" % (self.host, self.port, len(self.threads))
    
    def open_socket(self):
        """Opens a TCP/IP socket. 
        
        Also a currently running socket on the same
        port is recycled. This is a workaround for the typical 'Address 
        already in use' bug"""
        print "I am opening the server socket ..."
        try:
            self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            # self.server.setblocking(0) 
            # TODO: Learn the difference of blocking and non-blocking sockets
            self.server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self.server.bind((self.host,self.port))
            self.server.listen(5)
            print "The server socket is open ..."
        except socket.error, (value,message):
            if self.server:
                self.server.close()
            print "Could not open socket: " + message
            sys.exit(1)
    
    def _killStopped(self, list):
        """Receives a list of threads and joins all the non-active ones. 
        
        It returns a list of the active threads."""
        t = []
        for i in list:
            if i.isAlive():
                t.append(i)
            else:
                i.join()
        return t
    
    def run(self):
        """Mainloop of the server."""
        self.open_socket()
        input = [self.server,sys.stdin]
        running = 1
        print "I am starting the server mainloop ... waiting for connections"
        while running:
            inputready,outputready,exceptready = select.select(input,[],[])
            for s in inputready:
                if s == self.server:
                    (self.sock, address) = self.server.accept()
                    # self.transaction()
                    c = ConnectionThread(Transaction(sock=self.sock,
                    outputQueue=self.outputQueue, inputQueue=self.inputQueue))
                    c.start()
                    self.threads.append(c)
                    self.threads = self._killStopped(self.threads)
                elif s == sys.stdin:
                    junk = sys.stdin.readline()
                    running = 0
        self.server.close()
        for c in self.threads:
            c.join()
        print "Server shutting down. Bye Bye."
        sys.exit()


class ConnectionThread(threading.Thread):
    """A simple thread class that handles transactions."""
    def __init__(self, transaction):
        super(ConnectionThread, self).__init__()
        self.transaction = transaction
    def run(self):
        """Starts the transaction."""
        print "Greetings from a thread. My Names is", self.getName()
        self.transaction.run()


class Client(Networking):
    """An abstract client class.
    
    Gets initialized with a preFunc pointer. The client gets the arguments
    for the function from the server. It returns the return value of the 
    function back to the server."""
    def __init__(self, task=None, *args, **kwargs):
        super(Client, self).__init__(*args, **kwargs)
        self.task = task
    
    def run(self):
        """docstring for run"""
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # self.sock.setblocking(0)
        try:
            self.sock.connect((self.host, self.port))
        except socket.error, e:
            print "cannot establish connection to server", e
            self.sock.close()
        try:
            # stuff = self.recv_end()
            stuff = self.recv_pickle()
        except socket.error, e:
            print e
            self.sock.close()
            sys.exit()
        else:
            print stuff
            stuff2 = self.task.run(stuff)
            # self.send_end(stuff2)
            # self.sock.sendall(stuff2)
            self.send_pickle(stuff2)            
            self.sock.shutdown(socket.SHUT_WR)
            self.sock.close()


if __name__=='__main__':
    print Networking(len).port


# class DataClient(threading.Thread):
#     def __init__(self, queue):
#         super(DataClient, self).__init__()
#         self.queue = queue
#         self.d = dict()
#         self.activityDict = dict()
#
#     def run(self):
#         counter = 0
#         when2write = range(1,201000,100)
#         while True:
#             data = self.queue.get()
#             if data == None:
#                 break
#             counter += 1
#             self.d = sumDicts(self.d, data[0])
#             self.activityDict = sumDicts(self.activityDict, data[1])
#             print "Hi I am the DataClient. My Name is", self.getName(), filterDict(self.d, 0.9)
#             if counter in when2write:
#                 prefix = "/home/niko/SimulationData/relativeEssentiality/"
#                 path = prefix + "dict" + str(counter) + ".txt"
#                 print "writing to " + path
#                 mathDict1 = dict2mathematica(self.d)
#                 mathDict2 = dict2mathematica(self.activityDict)
#                 open(path, 'w').write(repr((mathDict1, mathDict2)))
