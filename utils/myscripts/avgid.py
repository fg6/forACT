from __future__ import print_function 
import sys
from numpy import genfromtxt
from numpy import mean

myfile=sys.argv[1]

chunkids,bases =genfromtxt(myfile,delimiter='\t', dtype=float, unpack=True, usecols =(2,3))
avgid=mean(chunkids)
cov=sum(bases)
print("%.2f" % avgid, "%d" % cov)

