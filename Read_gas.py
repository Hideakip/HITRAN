__author__ = 'Hideaki'
import numpy as np
import time

filename="C:\\Users\\Hideaki\\Dropbox\\Scytronix\\Software\\Eclipse projects\\HITRAN v2\\Raw data par\\CO2.npy"
t0=time.time()
a=np.load(filename)
t1=time.time()-t0
print "File read ",t1
#b=np.argsort(a,axis=0)

d=a[:,0]
t2=time.time()-t0
print "Sort=",t2-t1
print np.size(a[:,1])
print d[0],d[-1]