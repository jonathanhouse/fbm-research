from data_files import DataFile 
import matplotlib.pyplot as plt 
import numpy as np
from plot import plot, msd_fit


p = DataFile("data/gamma=1.0/nt=2**23/weight=-2.0/L=1000/intv=[0,2**23]")
p1 = DataFile("data/gamma=1.0/nt=2**23/weight=-2.0/L=1000/intv=[100]")

print(p.path,p.dis['P(x)*L'],p.length)
print(p.avx['<r^2>'])

fig, ax = plt.subplots(1,2)

plot(ax=ax, data=[[p],[p1]], type=['dis','msd'], label=['t','weight'])
msd_fit(ax=ax[1], data=[p1], interval=[96,193])

plt.show()