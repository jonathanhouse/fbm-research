from parse_data import get_data,get_tempered_data, get_err_data, get_fast_data
import matplotlib.pyplot as plt 
import numpy as np
import tikzplotlib as tz
from scipy.interpolate import CubicSpline
from numpy.polynomial.polynomial import Polynomial as Poly


w1_avx, w1_dis = get_data("data/gamma=1.0/nt=2**23/weight=-1.0/L=10M")
w0_avx, w0_dis = get_data("data/gamma=1.0/nt=2**23/weight=0.0/L=10M")
w025_avx, w025_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.25/L=10M")

w5_avx, w5_dis = get_data("data/gamma=1.0/nt=2**23/weight=-5.0/L=10M")


fig,ax = plt.subplots()
#x,y = 'x/L','P(x)*L' 
label1,label2,label3, label4,label5 = 'none', 'orig','both','$\gamma$ = 0.7','$\gamma$ = 0.2'

title = "different memory dis at $\gamma=0.2$" # title of plot 
sub = 'NT = 16384' # subtitle of plot 
#ax[0:2].set(xlabel='$'+x+'$',ylabel='$'+y+'$') # sets axis labels 
#ax[0:2].set(xscale='log',yscale='log') # sets axis scales 
#ax[0:2].set(title=title + '\n' + sub) # sets title/subtitle 

ax.set(xlabel='$t$',ylabel='$ \\langle x \\rangle $',title='free-space msd w/ power law\n$\gamma=1.0$, nt=2**23, L=10M, nconf=20000') # sets axis labels 

ax.set(xscale='log',yscale='log')

# power law stuff 
'''
x,y = np.log10(w1_avx["t"]), np.log10(w1_avx["<r^2>"])
print(x)
low_bound = np.nonzero(np.fabs(x - 1) < 1e-10)[0][0]
up_bound = np.nonzero(np.fabs(x - 2.05690485) < 1e-5)[0][0]
series = poly.Polynomial.fit(x[low_bound:up_bound],y[low_bound:up_bound],deg=1)

series = series.convert()
print(series.coef)

ax.plot(x,series(x))
ax.plot(x,y)

x_power = np.linspace(1,2,100)
ax.plot(10**x_power, 10**series(x_power))
'''

bound_arr = w1_avx['t']
x,y = np.log10(w1_avx['t']), np.log10(w1_avx['<r^2>'])
low_bound = np.nonzero(np.fabs( bound_arr - 3775787) < 1e-10)[0][0]
high_bound = np.nonzero(np.fabs( bound_arr - 7551573) < 1e-10)[0][0]

p1 = Poly.fit(np.log10(w1_avx['t'][low_bound:high_bound]), np.log10(w1_avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
p025 = Poly.fit(np.log10(w025_avx['t'][low_bound:high_bound]), np.log10(w025_avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
p5 = Poly.fit(np.log10(w5_avx['t'][low_bound:high_bound]), np.log10(w5_avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
p0 = Poly.fit(np.log10(w0_avx['t'][low_bound:high_bound]), np.log10(w0_avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)

p1 = p1.convert()
p0 = p0.convert()
p025 = p025.convert()
p5 = p5.convert()

x_test = np.linspace(np.log10(3775787), np.log10(7551573), 100)

ax.plot(w0_avx['t'],w0_avx['<r^2>'],label='weight=0.0\n' + str(p0.coef[1]))
ax.plot(w025_avx['t'],w025_avx['<r^2>'],label='weight=-0.25 \n' + str(p025.coef[1]))
ax.plot(w1_avx['t'],w1_avx['<r^2>'],label='weight=-1.0\n' + str(p1.coef[1]))
ax.plot(w5_avx['t'],w5_avx['<r^2>'],label='weight=-5.0\n' + str(p5.coef[1]))

ax.plot(10**x_test, 10**p1(x_test))
ax.plot(10**x_test, 10**p025(x_test))
ax.plot(10**x_test, 10**p5(x_test))
ax.plot(10**x_test, 10**p0(x_test))
#ax.plot(w1_avx['t'],w1_avx['<r^2>'],label='weight=-1.0')




#ax.plot(w016_dis['x/L'],w016_dis['P(x)*L'],label='weight=0.16')
#ax.plot(w5_dis['x/L'],w5_dis['P(x)*L'],label='weight=5.0')
#ax.plot(w1_dis['x/L'],w1_dis['P(x)*L'],label='weight=0.0')
#ax.plot(w016_avx[x],w016_avx[y],label='weight=-0.16')
#ax.plot(w1_avx[x],w1_avx[y],label='weight=-1.0')

#ax.plot(w5_avx[x],w5_avx[y],label='weight=-5.0')

'''
ax.plot(w016_avx[x],w016_avx[y],label='weight=-0.16')
#x.plot(w05_avx[x],w05_avx[y],label='weight=-0.5')
#ax.plot(w2_avx[x],w2_avx[y],label='weight=-2.0')

ax.plot(w5_avx[x],w5_avx[y],label='weight=-5.0')
'''

ax.legend()

plt.savefig("plot.svg")
plt.show()
