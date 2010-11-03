from tables import *
from numpy import *
h5 = openFile('test3.h5')
b=h5.root.simulations[:]
b['lowerBounds'][0]
set(b['lowerBounds'][0])
