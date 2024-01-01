from parse_data import get_data,get_tempered_data, get_err_data, get_fast_data
import matplotlib.pyplot as plt 
import numpy as np
from data_files import DataFile
from numpy.polynomial.polynomial import Polynomial as Poly


def plot_binned(ax,fig,data, binsize,label,type='linear',linestyle=None):
    L = data[0].length
    for d in data: 

       # if(d.length != L):
        #    d.dis["P(x)"] = np.array( list(np.zeros(shape=int((L - d.length)/2))) + list(d.dis["P(x)"]) + list(np.zeros(shape=int((L - d.length)/2))) )
        #    d.nbin = data[0].nbin
         #   d.dis["ibin"] = data[0].dis["ibin"]

        binned, bin_mid = ordered_binning(df=d,binsize=binsize)
        return binned,bin_mid
        if type == 'linear':
            ax.set(xscale='linear',yscale='linear')
            ax.set(xlabel='$x$',ylabel='$P(x)$')
        if type == 'log':
            ax.set(xscale='linear',yscale='log')
            ax.set(xlabel='$x$',ylabel='$P(x)$')
        if type == 'gaus':
            ax.set(xscale='linear',yscale='log')
            ax.set(xlabel='$x^2$',ylabel='$P(x)$')
            bin_mid = bin_mid**2*np.sign(bin_mid)
        ax.plot(bin_mid,binned,label=d.get_label(label),linestyle=linestyle)
    ax.legend()

    cols = len(data)
    title = ["nconf","length","gamma","nt","weight","t"]
    sup = "bin width=" + str(int(binsize)) + ", "
    for t in title:
        if t not in label: 
            sup += str(data[cols-1].get_label(t)) + ", "
    sup = sup[:-2]
    fig.suptitle(t=sup)



def ordered_binning(df, binsize):
    px = df.dis["P(x)"]     
    partitions = int(df.nbin*2/binsize)
    binned_dis = np.zeros(shape=int(partitions))
    midpoint_bins = np.zeros(shape=int(partitions))
    for n in range(int(partitions)):
        bin_start = binsize*n
        binned_dis[n] = np.sum( px[bin_start:bin_start + binsize] )
        midpoint_bins[n] = df.dis["ibin"][bin_start + int(0.5*binsize)]

    return binned_dis, midpoint_bins

def plot_times(ax,times,top):
    for t in times:
        ax.plot(t*np.ones(10),np.linspace(0,top,10))

def msd_fit(ax, data, interval):

    bound_arr = data[0].avx['t']
    low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-10)[0][0]
    high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-10)[0][0]
    for p in data: 



        x_test = np.linspace(np.log10(interval[0]),np.log10(interval[1]),100)

        series = Poly.fit(np.log10(p.avx['t'][low_bound:high_bound]), np.log10(p.avx['<r^2>'][low_bound:high_bound]), deg=1, window=None)
        series = series.convert()

        p.series = series

        ax.plot(10**x_test, 10**series(x_test),marker='x',markersize=3)
        ax.legend()


def log_fit(ax,data_x,data_y,interval):
    bound_arr = data_x
    low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-7)[0][0]
    high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-7)[0][0]

    x_test = np.linspace(np.log10(interval[0]),np.log10(interval[1]),100)
    series = Poly.fit(np.log10(data_x[low_bound:high_bound]), np.log10(data_y[low_bound:high_bound]), deg=1, window=None)
    series = series.convert()

    ax.plot(10**x_test, 10**series(x_test),marker='x',markersize=3,label=series)
    ax.legend()

    return series 



def plot(ax, fig, data, type, label,marker=None,xlim=None,markersize=None,linestyle=None,ylim=None):

    if(isinstance(ax,plt.Axes)): 
       
        cols = len(data)
        if(type == 'msd'): 
            ax.set(xscale = 'log', yscale = 'log')
            for p in data: 

                s = ""
                for l in label:
                    s += p.get_label(l) + ", "
                s = s[:-2]

                ax.plot(p.avx['t'], p.avx['<r^2>'],label=s,marker=marker,linestyle=linestyle)
                ax.set(xlabel='$t$',ylabel='$\\langle x^2 \\rangle$')

        if (type == 'dis'):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            for p in data: 
                
                s = ""
                for l in label:
                    s += p.get_label(l) + ", "
                s = s[:-2]

                ax.plot(p.dis['x'], p.dis['P(x)'], label=s,marker=marker,markersize=markersize,linestyle=linestyle)
                ax.set( xlabel='$x$',ylabel='$P(x)$')


        if (type == 'gaus'):
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set(xscale = 'linear', yscale = 'log')
            for p in data: 
                
                s = ""
                for l in label:
                    s += p.get_label(l) + ", "
                s = s[:-2]

                ax.plot(p.dis['x']**2*np.sign(p.dis['x']), p.dis['P(x)'], label=s,marker=marker,markersize=markersize,linestyle=linestyle)
                ax.set( xlabel='$x^{2}$',ylabel='$P(x)$')


        ax.legend()

        title = ["nconf","length","gamma","nt","weight","t"]
        sup = ""
        for t in title:
            if t not in label: 
                sup += str(data[cols-1].get_label(t)) + ", "
        sup = sup[:-2]
        fig.suptitle(t=sup)


    else:  

        cols = len(ax)
        for i in range(cols):

            if(type[i] == 'msd'): 
                ax[i].set(xscale = 'log', yscale = 'log')
                for p in data[i]: 
                    ax[i].plot(p.avx['t'], p.avx['<r^2>'],label= p.get_label(label[i]))
                    ax[i].set(xlabel='$t$',ylabel='$\\langle x^2 \\rangle$')

            if (type[i] == 'dis'):
                ax[i].set_ylim(ylim)
                ax[i].set_xlim(xlim)
                for p in data[i]: 
                    ax[i].plot(p.dis['x'], p.dis['P(x)'], label=p.get_label(label[i]))
                    ax[i].set( xlabel='$x$',ylabel='$P(x)$')

            if (type[i] == 'gaus'):
                ax[i].set_xlim(xlim)
                ax[i].set_ylim(ylim)
                ax[i].set(xscale = 'linear', yscale = 'log')
                for p in data[i]: 

                    ax[i].plot(p.dis['x']**2*np.sign(p.dis['x']), p.dis['P(x)'], label=p.get_label(label[i]))
                    ax[i].set( xlabel='$x^{2}$',ylabel='$P(x)$')

            ax[i].legend()

        title = ["nconf","length","gamma","nt","weight","t"]
        sup = ""
        for t in title:
            if t not in label: 
                sup += str(data[cols-1][0].get_label(t)) + ", "
        sup = sup[:-2]
        fig.suptitle(t=sup)