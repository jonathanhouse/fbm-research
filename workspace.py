from data_files import DataFile 
import matplotlib.pyplot as plt 
import numpy as np
from plot import plot, msd_fit, plot_times, plot_binned
from numpy.polynomial.polynomial import Polynomial as Poly
import tikzplotlib as tz 


'''
x1 = DataFile("data/linear force/gamma=1.0/weight=-unifom[0,0.5]/nt=2**26/L=100/asymmetric gradient/")
x2 = DataFile("data/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=100/mixed gradient/nbin=50/")
x3 = DataFile("data/nonlinear force/prefactor=0.25/gamma=1.0/weight=-0.25/L=100/asymmetric gradient/nbin=100/")
x4 = DataFile("data/nonlinear force/prefactor=1.0/gamma=1.0/weight=-0.25/nt=2**26/L=100/nbin=50/mixed gradient/")
x5 = DataFile("data/linear force/gamma=1.0/weight=0.0/nt=2**23/L=10M")

fig,ax =  plt.subplots(1,1)

ax.plot(x1.dis["x/L"],x1.dis["P(x)*L"],label='-uniform[0,0.5],asym,nbin=50')
ax.plot(x2.dis["x/L"],x2.dis["P(x)*L"],label='weight=-0.25,mixed gradient,nbin=50')
ax.plot(x3.dis["x/L"],x3.dis["P(x)*L"],label='nonlinear=0.25,weight=-0.25,asym,nbin=100')
ax.plot(x4.dis["x/L"],x4.dis["P(x)*L"],label='nonlinear=1.0,weight=-0.25,mixed,nbin=50')

ax.set(xlabel='x/L',ylabel='P(x)*L')


fig.suptitle("nconf=25K, L=100, $\gamma=1.0$")

ax.legend()
'''

x10 = DataFile("../bfbm/data/organized_bfbm/orig/gamma=1.3",type='grad')

fig,ax = plt.subplots(1,1)

interval=[0.106280E+01, 0.121193E+01]
bound_arr = x10.log['x']
low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-10)[0][0]
high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-10)[0][0]
p = x10

x_test = np.linspace(np.log10(interval[0]),np.log10(interval[1]),100)

series = Poly.fit(np.log10(p.log['x'][low_bound:high_bound]), np.log10(p.log['P(x)'][low_bound:high_bound]), deg=1, window=None)
print(series)
series = series.convert()
print(series)
p.series = series

ax.plot(10**x_test, 10**series(x_test),marker='x',markersize=3)



ax.plot(x10.log["x"],x10.log["P(x)"],label=x10.series)
x =np.linspace(0.213052E+01,0.590629E+01,100)
y = -0.78932762*x_test-(1.84380828)


y2 = -2.26197158 - 0.15888054*x
#ax.plot(x,10**y,'ro')
#ax.plot(x,10**y2,'bo')
#ax.plot(x,(x**( -0.78932762))*(0.0143282028051),'go')

ax.set(xscale='log',yscale='log')
ax.legend()


#ax.plot(x10.dis["x/L"],x10.dis["P(x)*L"],label='nbin=' + str(x10.nbin))
#ax.plot(x20.dis["x/L"],x20.dis["P(x)*L"], label= 'nbin=' + str(x20.nbin))
#ax.set(xlim=[-0.5,0.5])


#ax.legend()

#plot(ax=ax,fig=fig,data=[x10],type='dis',label=['nbin'])

#ax.set(title='comparing symmetric linear fit gradient: window=3')

#ax.legend()
'''

x = DataFile("data/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M/linear fit gradient/symmetric gradient/window=3/nbin=5M")
x2 = DataFile("data/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M/linear fit gradient/symmetric gradient/window=3/nbin=10M")
x3 = DataFile("data/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M/linear fit gradient/symmetric gradient/window=3/nbin=20M")
fig,ax = plt.subplots(1,1)

ax.plot(x.avx['t'],x.avx['<r^2>'],label='nbin=5M')
ax.plot(x2.avx['t'],x2.avx['<r^2>'],label='nbin=10M')
ax.plot(x3.avx['t'],x3.avx['<r^2>'],label='nbin=20M')

ax.set(xscale='log',yscale='log')


fig.suptitle('nt='+str(x.nt) + ',nconf='+str(x.nconf)+ ',weight='+str(x.weight) + ",L="+str(x.length) + ",$\gamma$=" +str(x.gamma))
ax.legend()

plt.show()
'''


'''
g13 = [ DataFile("../bfbm/data/22-23 data/02:22/orig/gamma=1.6",type='grad'),
       DataFile("../bfbm/data/22-23 data/02:22/none/gamma=1.6",type='grad'),
       DataFile("../bfbm/data/organized_bfbm/both/gamma=1.6",type='grad'),
]

g02 = [ DataFile("../bfbm/data/22-23 data/02:22/orig/gamma=0.2",type='grad'),
       DataFile("../bfbm/data/22-23 data/02:22/none/gamma=0.2",type='grad'),
       DataFile("../bfbm/data/22-23 data/02:22/both/gamma=0.2",type='grad'),
]
'''
'''
g13 = [ DataFile("../bfbm/data/22-23 data/02:22/orig/gamma=1.6",type='grad'),
       DataFile("../bfbm/data/organized_bfbm/none/gamma=1.6",type='grad'),
       DataFile("../bfbm/data/organized_bfbm/both/gamma=1.6",type='grad'),
]

g02 = [ DataFile("../bfbm/data/fall-23-runs/orig/gamma=0.2",type='fast hosking'),
       DataFile("../bfbm/data/fall-23-runs/none//gamma=0.2",type='fast hosking'),
       DataFile("../bfbm/data/22-23 data/02:22/both/gamma=0.2",type='grad'),
]

fig,ax = plt.subplots(2,1)
ax[1].set(xscale='log',yscale='log')
ax[0].set(xscale='log',yscale='log')

ax[0].set(xlabel='$x/L$',ylabel='$P(x)L$',title='$\\alpha=1.8$')
ax[1].set(xlabel='$x/L$',ylabel='$P(x)L$',title='$\\alpha=0.4$')

ax[1].plot(g13[0].log["x/L"],g13[0].log["P(x)*L"],label='ORIG')
ax[1].plot(g13[1].log["x/L"],g13[1].log["P(x)*L"],label='NONE')
ax[1].plot(g13[2].log["x/L"],g13[2].log["P(x)*L"],label='BOTH')

ax[0].plot(g02[0].log["x/L"],g02[0].log["P(x)*L"],label='ORIG')
ax[0].plot(g02[1].log["x/L"],g02[1].log["P(x)*L"],label='NONE')
ax[0].plot(g02[2].log["x/L"],g02[2].log["P(x)*L"],label='BOTH')

ax[0].legend()
ax[1].legend()
'''

plt.show()






#plot(ax=ax,fig=fig,type=['dis','dis'],label=['gamma','gamma'], data=[g02,g13])

'''

print(np.shape(d))
print(np.shape(x1))
print(np.shape(x2))

A = np.array([x1,x2]).T

print(np.shape(A))

pA = np.linalg.pinv(A)

print( pA )
print(pA @ d.T)

x = pA @ d.T

x_plot = np.arange(-5e6,5e6+1,1)
y_plot = A @ x

plot(ax=ax,fig=fig,type='dis',label=["weight"],data=cumdis)
ax.plot(x_plot,y_plot,linestyle="dashed")
ax.set_xlim([-1e6,1e6])
ax.set_ylim([0,1e-4])
'''


#plt.show()
