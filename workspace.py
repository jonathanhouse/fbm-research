from data_files import DataFile 
import matplotlib.pyplot as plt 
import numpy as np
from plot import plot, msd_fit, plot_times, plot_binned, log_fit
from numpy.polynomial.polynomial import Polynomial as Poly
import tikzplotlib as tz 


fig,ax = plt.subplots(1,1)
y = DataFile("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**10/L=10M (nbin=5M)/nconf=64*32/asymmetric gradient",'avx')
x = DataFile("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**10/L=10M (nbin=5M)/nconf=64*512/asymmetric gradient",'avx')

x1 = DataFile("data/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=1.5M/intv=[2**26]",'avx')

ax.set(xscale='log',yscale='log')
ax.set(xlabel='x/L')


fig.suptitle('nbin/L=5M/10M, weight=-0.25, $\gamma=1.0$, asymmetric gradient')
#ax.plot(x.cor['t'],(x.cor['<r^2>']))

pos = y.avx
pos1 = x.avx
pos3 = x1.avx


ax.plot(pos['t'],pos['<r^2>'],label='nconf=64*32')
#ax.plot(pos1['t'],pos1['<r^2>'],label='nconf=64*512')

ax.plot(pos3['t'],pos3['<r^2>'],label='nconf=10*64,non-parallel')

#log_fit(ax,pos['t'],pos['<r^2>'],interval=[461,922])


ax.legend()
#log_fit(ax,g1_w5_steps.cor['t'],g1_w5_steps.cor['<|grad_pos|>'],interval=[8980384,60412582])
'''


ax.plot(y.cor['t'],(y.cor['<r^2>']),label='asymmetric gradient <r^2>',linestyle='--' )
ax.plot(y.cor['t'],(y.cor['<xix_pos^2>']),label='asymmetric gradient <xix_pos^2>' ,linestyle='--')
ax.plot(y.cor['t'],(y.cor['<f_grad_pos^2>']),label='asymmetric gradient <f_grad_pos^2>' ,linestyle='--')
ax.plot(y.cor['t'],(y.cor['<xix_pos*f_grad_pos>']),label='asymmetric gradient <xix_pos*f_grad_pos>' ,linestyle='--')
'''

'''
ax.plot(x.cor['t'],(x.cor['<r^2>']),label='<r^2>')
ax.plot(x.cor['t'],(x.cor['<xix_pos^2>']),label='<xix_pos^2>')
ax.plot(x.cor['t'],(x.cor['<f_grad_pos^2>']),label='<f_grad_pos^2>' )
ax.plot(x.cor['t'],(x.cor['<xix_pos*f_grad_pos>']),label='<xix_pos*f_grad_pos>')

ax.plot(x1.cor['t'],(x1.cor['<|xix_pos|>']),label='<|xix_step|>' )
ax.plot(x1.cor['t'],(x1.cor['<|grad_pos|>']),label='<|f_grad_step|>' )

y = np.linspace(1e-3,1e11,10)
ax.plot(22878*np.ones(np.size(y)),y,label='t=22878; <|xix_step|> & <|f_grad_step|> intersection')

y = np.linspace(1e-3,1e11,10)
ax.plot(6.75e6*np.ones(np.size(y)),y,label='t=6.75e6; <xix_pos^2> & <f_grad_pos^2> intersection')
'''

#series = log_fit(ax,x1.cor['t'],x1.cor['<|grad_pos|>'],interval=[461,41717])
#series = log_fit(ax,x1.cor['t'],x1.cor['<|grad_pos|>'],interval=[6350091,60412582])

'''
ax.plot(y.cor['t'],(y.cor['<xix_pos^2>']),label='<xix_pos^2>; symmetric')
ax.plot(y.cor['t'],(y.cor['<f_grad_pos^2>']),label='<f_grad_pos^2>; symmetric')
ax.plot(y.cor['t'],(y.cor['<xix_pos*f_grad_pos>']),label='<xix_pos*f_grad_pos>; symmetric')

ax.plot(y1.cor['t'],(y1.cor['<|xix_pos|>']),label='<|xix_pos|>; symmetric' )
ax.plot(y1.cor['t'],(y1.cor['<|grad_pos|>']),label='<|f_grad_pos|>; symmetric' )

'''




'''



L100M_data = DataFile("data/linear force/gamma=0.4/weight=-0.25/nt=2**26/L=100M/nbin=50M/asymmetric gradient/intv=[2**26]",'avx')
#L10M_asym_data = DataFile("data/linear force/gamma=0.4/weight=-0.25/nt=2**26/L=10M/")
L10M_sym_data = DataFile("data/linear force/gamma=0.4/weight=-0.25/nt=2**26/L=10M/nbin=5M/asymmetric gradient",'avx')


bm_data = DataFile("data/linear force/gamma=0.4/weight=0.0/nt=2**26/L=10M")



fig,ax = plt.subplots(1,1)

fig.suptitle("nt=" + str(bm_data.nt) + ", nconf=" + str(bm_data.nconf) + ", $\gamma$=0.4" + ", weight=" + str(L100M_data.weight))

msd_fit(ax,[L100M_data],[25400363,42718147])
msd_fit(ax,[bm_data],[25400363,42718147])

ax.set(xscale='log',yscale='log')

ax.plot(L100M_data.avx["t"],L100M_data.avx["<r^2>"],label='L=100M,asymmetric gradient : ' + str(L100M_data.series))
#ax.plot(L10M_asym_data.avx["t"],L10M_asym_data.avx["<r^2>"],label='L=10M,asymmetric gradient')
ax.plot(L10M_sym_data.avx["t"],L10M_sym_data.avx["<r^2>"],label='L=10M,asymmetric gradient')
ax.plot(bm_data.avx["t"], bm_data.avx["<r^2>"],label='pure FBM : ' + str(bm_data.series)) 

ax.set(xlabel='t',ylabel='<r^2>')
'''

ax.legend()
'''
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
'''


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
