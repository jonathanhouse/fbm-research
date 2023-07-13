from parse_data import get_data,get_tempered_data, get_err_data, get_fast_data
import matplotlib.pyplot as plt 
import numpy as np
from data_files import DataFile
from numpy.polynomial.polynomial import Polynomial as Poly

def plot_times(ax,times,top):
    for t in times:
        ax.plot(t*np.ones(10),np.linspace(0,top,10))

def msd_fit(ax, data, interval):

    bound_arr = data[0].avx['t']
    for p in data: 

        low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-10)[0][0]
        high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-10)[0][0]
        x_test = np.linspace(np.log10(interval[0]),np.log10(interval[1]),100)

        series = Poly.fit(np.log10(p.avx['t'][low_bound:high_bound]), np.log10(p.avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
        series = series.convert()

        ax.plot(10**x_test, 10**series(x_test),label=str(series))
        ax.legend()


def plot(ax, fig, data, type, label):

    if(isinstance(ax,plt.Axes)): 

        if(type == 'msd'): 
            ax.set(xscale = 'log', yscale = 'log')
            for p in data: 
                ax.plot(p.avx['t'], p.avx['<r^2>'],label=p.get_label(label))
                ax.set(xlabel='$t$',ylabel='$\\langle x^2 \\rangle$')

        if (type == 'dis'):
            for p in data: 
                ax.plot(p.dis['x/L'], p.dis['P(x)*L'], label= p.get_label(label))
                ax.set( xlabel='$x/L$',ylabel='$P(x) \\times L$')

        ax.legend()


    else:  

        cols = len(ax)
        for i in range(cols):

            if(type[i] == 'msd'): 
                ax[i].set(xscale = 'log', yscale = 'log')
                for p in data[i]: 
                    ax[i].plot(p.avx['t'], p.avx['<r^2>'],label= p.get_label(label[i]))
                    ax[i].set(xlabel='$t$',ylabel='$\\langle x^2 \\rangle$')

            if (type[i] == 'dis'):
                for p in data[i]: 
                    ax[i].plot(p.dis['x/L'], p.dis['P(x)*L'], label=p.get_label(label[i]))
                    ax[i].set( xlabel='$x/L$',ylabel='$P(x) \\times L$')

            ax[i].legend()

        title = ["nconf","length","gamma","nt","weight","t"]
        sup = ""
        for t in title:
            if t not in label: 
                sup += str(data[cols-1][0].get_label(t)) + ", "
        sup = sup[:-2]
        fig.suptitle(t=sup)