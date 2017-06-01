import sys
from numpy import genfromtxt

#def split(x):
#	return x.split("_")[0]


myfile=sys.argv[1]

l = genfromtxt(myfile,delimiter='\t', dtype=str, unpack=True, usecols =0)

#intl=[int(x) for x in map(split, l)]

l1=(list(set(l)))

for ctg in l1:
	print(ctg)
