import subprocess as sub
from numpy import *


runs = 100

xs = []

for _ in range(runs):
    proc = sub.Popen("level-n-trees.exe", stdout = sub.PIPE)
    out = int(proc.communicate()[0])
    print out
    xs.append(out)
    
print xs
print mean(xs)
print var(xs)
    