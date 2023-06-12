from parse_data import get_data,get_tempered_data, get_err_data, get_fast_data
import matplotlib.pyplot as plt 
import numpy as np
import tikzplotlib as tz
from scipy.interpolate import CubicSpline


w1_avx, w1_dis = get_data("data/gamma=0.4/nt=2**23/weight=-1.0/L=100000")
w0_avx, w0_dis = get_data("data/gamma=0.4/nt=2**23/weight=0.0/L=100000")
'''

w05_avx, w05_dis = get_data("data/gamma=0.4/nt=2**23/weight=-0.5/L=100000")
w2_avx, w2_dis = get_data("data/gamma=0.4/nt=2**23/weight=-2.0/L=100000")
w4_avx, w4_dis = get_data("data/gamma=0.4/nt=2**23/weight=-4.0/L=100000")
w25_avx, w25_dis = get_data("data/gamma=0.4/nt=2**23/weight=-0.25/L=100000")
'''

fig,ax = plt.subplots()
#x,y = 'x/L','P(x)*L' 
x,y = 't','<r^2>'
label1,label2,label3, label4,label5 = 'none', 'orig','both','$\gamma$ = 0.7','$\gamma$ = 0.2'

title = "different memory dis at $\gamma=0.2$" # title of plot 
sub = 'NT = 16384' # subtitle of plot 
#ax[0:2].set(xlabel='$'+x+'$',ylabel='$'+y+'$') # sets axis labels 
#ax[0:2].set(xscale='log',yscale='log') # sets axis scales 
#ax[0:2].set(title=title + '\n' + sub) # sets title/subtitle 

ax.set(xlabel='$t$',ylabel='$ \langle x^2 \\rangle$',title='$\gamma=0.4$ \nnt=2**23, L=100000, nconf=50000') # sets axis labels 

ax.set(xscale='log',yscale='log')

ax.plot(w0_avx[x],w0_avx[y],label='weight=0')
ax.plot(w1_avx[x],w1_avx[y],label='weight=-1')

'''

ax.plot(w25_avx[x],w25_avx[y],label='weight=-0.25')
ax.plot(w05_avx[x],w05_avx[y],label='weight=-0.5')

ax.plot(w2_avx[x],w2_avx[y],label='weight=-2')
ax.plot(w4_avx[x],w4_avx[y],label='weight=-4')
'''

ax.legend()

plt.savefig("plot.svg")
plt.show()
