from parse_data import get_data,get_tempered_data, get_err_data, get_fast_data
import matplotlib.pyplot as plt 
import numpy as np
import tikzplotlib as tz
from scipy.interpolate import CubicSpline
from numpy.polynomial.polynomial import Polynomial as Poly


gm1_avx, gm1_dis = get_data("data/gamma=1.0/nt=2**26/weight=-0.25/L=10M")
gm16_avx, gm16_dis = get_data("data/gamma=1.6/nt=2**26/weight=-0.25/L=10M")
gm04_avx, gm04_dis = get_data("data/gamma=0.4/nt=2**26/weight=-0.25/L=10M")
gm07_avx, gm07_dis = get_data("data/gamma=0.7/nt=2**26/weight=-0.25/L=10M")
gm13_avx, gm13_dis = get_data("data/gamma=1.3/nt=2**23/weight=-1.0/L=1000/intv=[100]")
gm11_avx, gm11_dis = get_data("data/gamma=1.1/nt=2**23/weight=-1.0/L=1000/intv=[100]")
gm06_avx, gm06_dis = get_data("data/gamma=0.6/nt=2**26/weight=-0.25/L=10M")
gm05_avx, gm05_dis = get_data("data/gamma=0.5/nt=2**26/weight=-0.25/L=10M")
gm08_avx, gm08_dis = get_data("data/gamma=0.8/nt=2**26/weight=-0.25/L=10M")
gm09_avx, gm09_dis = get_data("data/gamma=0.9/nt=2**26/weight=-0.25/L=10M")

gen_avx, gen_dis = get_data("data/gamma=1.0/nt=2**23/weight=-2.0/L=1000/intv=[0,2**23]")

w0_avx, w0_dis = get_data("data/gamma=1.0/nt=2**23/weight=0.0/L=1000")
w1_avx, w1_dis = get_data("data/gamma=1.6/nt=2**26/weight=-1.0/L=10M")
w5_avx, w5_dis = get_data("data/gamma=0.4/nt=2**23/weight=-5.0/L=1000/intv=[100]")
w016_avx, w016_dis = get_data("data/gamma=0.4/nt=2**23/weight=-0.16/L=1000/intv=[100]")
w025_avx, w025_dis = get_data("data/gamma=1.6/nt=2**26/weight=-0.25/L=10M")
w005_avx, w005_dis = get_data("data/gamma=1.6/nt=2**26/weight=-0.05/L=10M")

intv100_avx, intv100_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,100]")
intv1k_avx, intv1k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,1000]")
intv7k_avx, intv7k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,7000]")
int1v4k_avx, intv14k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,14000]")
intv40k_avx, intv40k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,40000]")
intv70k_avx, intv70k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,70000]")
intv160k_avx, intv160k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[0,160000]")

t100_avx, t100_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[100]")
t1k_avx, t1k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[1000]")
t7k_avx, t7k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[7000]")
t14k_avx, t14k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[14000]")
t40k_avx, t40k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[40000]")
t70k_avx, t70k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[70000]")
t160k_avx, t160k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=1000/intv=[160000]")

L100k_avx, L100k_dis = get_data("data/gamma=1.0/nt=2**26/weight=-0.5/L=100K")

intv100k_avx, intv100k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=100K/intv=[0,10**5]")
t100k_avx, t100k_dis = get_data("data/gamma=1.0/nt=2**23/weight=-0.5/L=100K/intv=[10**5,10**5]")

fig,ax = plt.subplots(1,3)
#x,y = 'x/L','P(x)*L' 
label1,label2,label3, label4,label5 = 'none', 'orig','both','$\gamma$ = 0.7','$\gamma$ = 0.2'

title = "different memory dis at $\gamma=0.2$" # title of plot 
sub = 'NT = 16384' # subtitle of plot 
#ax[0:2].set(xlabel='$'+x+'$',ylabel='$'+y+'$') # sets axis labels 
#ax[0:2].set(xscale='log',yscale='log') # sets axis scales 
#ax[0:2].set(title=title + '\n' + sub) # sets title/subtitle 

fig.suptitle('$\gamma$=1.0, nt=2**23, L=100K, nconf=50000')

ax[2].set(xlabel='$t$',ylabel='$\\langle x^2 \\rangle$', title = "varying weights of msd")
ax[0].set(xlabel='$x/L$',ylabel='$P(x) \\times L$', title = 'instantaneous distributions of weight=-0.5')
ax[1].set(xlabel='$x/L$',ylabel='$P(x) \\times L$', title = 'cummulative distributions of weight=-0.5')


ax[0].set(xscale='linear',yscale='linear')
ax[1].set(xscale='linear',yscale='linear')
ax[2].set(xscale='log',yscale='log')

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



bound_arr = gm04_avx['t']

low_bound = np.nonzero(np.fabs( bound_arr - 30206292) < 1e-10)[0][0]
high_bound = np.nonzero(np.fabs( bound_arr - 42718147) < 1e-10)[0][0]

p04 = Poly.fit(np.log10(gm04_avx['t'][low_bound:high_bound]), np.log10(gm04_avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
p04 = p04.convert()
ax.plot(gm04_avx['t'], gm04_avx['<r^2>'],label='$\gamma$=0.4\n' + str(p04))
x_test = np.linspace(np.log10(30206292), np.log10(42718147), 100)
ax.plot(10**x_test, 10**p04(x_test))


'''

bound_arr = L100k_avx['t']

low_bound = np.nonzero(np.fabs( bound_arr - 70160) < 1e-10)[0][0]
high_bound = np.nonzero(np.fabs( bound_arr - 117994) < 1e-10)[0][0]

p04 = Poly.fit(np.log10(L100k_avx['t'][low_bound:high_bound]), np.log10(L100k_avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
p04 = p04.convert()

x_test = np.linspace(np.log10(70160), np.log10(117994), 100)



ax[2].plot(L100k_avx['t'],L100k_avx['<r^2>'],label='weight=-0.5\n' + str(p04))
#ax[2].plot(w0_avx['t'],w0_avx['<r^2>'],label='weight=0.0')

ax[2].plot(10**x_test, 10**p04(x_test))

ax[1].plot(intv100k_dis['x/L'],intv100k_dis['P(x)*L'],label='t=[0,10**5]')
ax[0].plot(t100k_dis['x/L'],t100k_dis['P(x)*L'],label='t=10**5')

ax[0].legend()
ax[1].legend()
ax[2].legend()

plt.savefig("plot.svg")
plt.show()
