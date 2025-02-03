from parse_data import get_data,get_tempered_data, get_err_data, get_fast_data
import matplotlib.pyplot as plt 
import numpy as np
from data_files import DataFile
from numpy.polynomial.polynomial import Polynomial as Poly
from decimal import Decimal
import scipy.optimize as opt

def suptitle_generator(x):
    title = "2*NBIN/L=" + "{:.2E}".format(Decimal(2*float(x.params["NBINS"]))) + "/" + "{:.2E}".format(Decimal(x.params["L"])) + "=" + str(2*float(x.params["NBINS"])/float(x.params["L"]))
    title += ", gamma=" + str(float(x.params["GAMMMA"])) 
    title += ", nconf=" + str(x.params["NCONF"]) 
    title += ", NT=$2^{" + str(np.log(float(x.params["NT"]))/np.log(2)) + "}$"
    return title






def calc_binsize(dis_in, num_bins):

    low_bound = 0
    high_bound = 0
    ep = 1e-10
    bin_per_x = 2*float(dis_in.params["NBINS"]) / float(dis_in.params["L"])
    high_bound = len(dis_in.dis["ibin"]) - 1
    for l in range(len(dis_in.dis["P(|x|)"])):
        if(dis_in.dis["P(|x|)"][l] > ep):
            low_bound = l
            break
    for h in range(len(dis_in.dis["P(|x|)"])):
        if(dis_in.dis["P(|x|)"][len(dis_in.dis["P(|x|)"]) - 1 - h] > ep):
            high_bound = len(dis_in.dis["P(|x|)"]) - 1 - h
            break
    binsize = int(bin_per_x * int((dis_in.dis["x"][high_bound] - dis_in.dis["x"][low_bound])/num_bins))
    print( "true binsize:", int(dis_in.dis["x"][high_bound] - dis_in.dis["x"][low_bound]) * bin_per_x / binsize )
    return binsize, low_bound, high_bound 

def unscaled_idis_plot(dis_list, scale_to_nt, rebinned=False, num_bins=1):
    fig,ax = plt.subplots(1,1)
    ax.set(xscale='linear',yscale='linear')
    ax.set(xlabel='x', ylabel='P(|x|)')
    fig.suptitle(suptitle_generator(dis_list[0]))
    title = suptitle_generator(dis_list[0])
    if rebinned: 
        title += ", num_bins=" + str(num_bins)
    title += "\n" + "unscaled idis"

    fig.suptitle(title)

    for i in range(len(dis_list)): 
        nt_ratio = scale_to_nt / float(dis_list[i].params["NT"]) 
        nt_range = float(dis_list[i].params["NTEND"]) - float(dis_list[i].params["NTSTART"])
        print("nt_range", nt_range)
        print("nt_ratio:", nt_ratio)
        # dis_list[i].dis["x"] = dis_list[i].dis["x"] * (nt_ratio)**(2.0/3.0)
        # dis_list[i].dis["P(|x|)"] = dis_list[i].dis["P(|x|)"] * (nt_ratio)**(-2.0/3.0)
        if(rebinned):
            bin_width, low, high = calc_binsize(dis_list[i], num_bins)
            print("bin_width:",bin_width)
            
            y_prime, bin_prime, x_prime = ordered_binning(df=dis_list[i],binsize=int(bin_width), low_bound=low, high_bound=high)
            y_prime = y_prime / nt_range / bin_width
            #y_prime = y_prime * (nt_ratio)**(-2.0/3.0)
            print("sum(y'):",sum(y_prime))
            print("len(y'):",len(y_prime))
            ax.plot(x_prime, y_prime, label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
        else:
            ax.plot(dis_list[i].dis["x"] , dis_list[i].dis["P(|x|)"], label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
    ax.legend()

def scaled_idis_plot(dis_list, scale_to_nt, rebinned=False, num_bins=1):
    fig,ax = plt.subplots(1,1)
    ax.set(xscale='linear',yscale='linear')
    ax.set(xlabel='x', ylabel='P(|x|)')
    fig.suptitle(suptitle_generator(dis_list[0]))
    title = suptitle_generator(dis_list[0])
    if rebinned: 
        title += ", num_bins=" + str(num_bins)
    title += "\n" + "scaled idis up to " + "$nt=2^{" + str(np.log(float(scale_to_nt))/np.log(2)) + "}$"

    fig.suptitle(title)

    for i in range(len(dis_list)): 
        nt_ratio = scale_to_nt / float(dis_list[i].params["NT"]) 
        nt_range = float(dis_list[i].params["NTEND"]) - float(dis_list[i].params["NTSTART"])
        print("nt_range", nt_range)
        print("nt_ratio:", nt_ratio)
        # dis_list[i].dis["x"] = dis_list[i].dis["x"] * (nt_ratio)**(2.0/3.0)
        # dis_list[i].dis["P(|x|)"] = dis_list[i].dis["P(|x|)"] * (nt_ratio)**(-2.0/3.0)
        if(rebinned):
            bin_width, low, high = calc_binsize(dis_list[i], num_bins)
            print("bin_width:",bin_width)
            
            y_prime, bin_prime, x_prime = ordered_binning(df=dis_list[i],binsize=int(bin_width), low_bound=low, high_bound=high)
            x_prime = x_prime * (nt_ratio)**(2.0/3.0)
            # y_prime = y_prime / nt_range 
            y_prime = y_prime * (nt_ratio)**(-2.0/3.0) / bin_width / nt_range
            print("sum(y'):",sum(y_prime) * (0.1 * bin_width) )
            print("len(y'):",len(y_prime))
            ax.plot(x_prime, y_prime, label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
        else:
            ax.plot(dis_list[i].dis["x"] , dis_list[i].dis["P(|x|)"], label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
    ax.legend()

def scaled_cdis_plot(dis_list, scale_to_nt, rebinned=False, num_bins=1):
    fig,ax = plt.subplots(1,1)
    ax.set(xscale='linear',yscale='linear')
    ax.set(xlabel='x', ylabel='P(|x|)')
    fig.suptitle(suptitle_generator(dis_list[0]))
    title = suptitle_generator(dis_list[0])
    if rebinned: 
        title += ", num bins under dis=" + str(num_bins)
    title += "\n" + "scaled cumulative distribution up to " + "$nt=2^{" + str(np.log(float(scale_to_nt))/np.log(2)) + "}$"
    fig.suptitle(title)

    markers = ['s', '*', 'o']
    for i in range(len(dis_list)): 
        nt_ratio = scale_to_nt / float(dis_list[i].params["NT"]) 
        # print(nt_ratio)

        if(rebinned):
            # dis_list[i].dis["P(|x|)"] = dis_list[i].dis["P(|x|)"] * (nt_ratio)**(1.0/3.0)
            bin_range, low_bound, high_bound  = calc_binsize(dis_list[i],num_bins)
            bin_range = 5e4 * (4/5.0)**i
            
            y_prime, bin_prime, x_prime = ordered_binning(df=dis_list[i],binsize=int(bin_range), low_bound=low_bound, high_bound=high_bound)
            x_prime = x_prime * (nt_ratio)**(2.0/3.0)
            y_prime = y_prime * (nt_ratio)**(1.0/3.0) / bin_range 
            print("sum(y')", sum(y_prime))
            print("len(y')", len(y_prime))
            ax.scatter(x_prime, y_prime, label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$", marker=markers[i])
        else:
            dis_list[i].dis["x"] = dis_list[i].dis["x"] * (nt_ratio)**(2.0/3.0)
            dis_list[i].dis["P(|x|)"] = dis_list[i].dis["P(|x|)"] * (nt_ratio)**(1.0/3.0)
            ax.plot(dis_list[i].dis["x"] , dis_list[i].dis["P(|x|)"], label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
    ax.legend()

def unscaled_dis_plot(dis_list, rebinned=False, num_bins=1):
    fig,ax = plt.subplots(1,1)
    ax.set(xscale='linear',yscale='linear')
    ax.set(xlabel='x', ylabel='P(|x|)')
    fig.suptitle(suptitle_generator(dis_list[0]))
    title = suptitle_generator(dis_list[0])
    if rebinned: 
        title += ", num bins under dis=" + str(num_bins)
    title += "\n" + "cumulative distributions time evolution"
    fig.suptitle(title)

    X = dis_list[0].dis["x"]
    low_bound = 0
    high_bound = 0
    ep = 1e-5
    for i in range(len(dis_list)): 
        if(rebinned):

            bin_range,low_bound, high_bound  = calc_binsize(dis_list[i],num_bins)
            
            y_prime, bin_prime, x_prime = ordered_binning(df=dis_list[i],binsize=int(bin_range),low_bound=low_bound, high_bound=high_bound )
            y_prime = y_prime / bin_range
            print("sum(y')", sum(y_prime))
            ax.plot(x_prime, y_prime, label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
        else:
            ax.plot(X, dis_list[i].dis["P(|x|)"], label="NT=$2^{" + str(np.log(float(dis_list[i].params["NT"]))/np.log(2)) + "}$")
            print("sum(P(|x|)')", sum(dis_list[i].dis["P(|x|)"]))
    ax.legend()



def signed_square(arr):
    return np.sign(arr)*np.power(arr,2)
def signed_sqrt(arr):
    return np.sign(arr)*np.power(np.abs(arr),0.5)


def plot_binned(ax,fig,data, binsize=1,label=None,type='linear',linestyle=None):
    #L = data[0].length

    for d in data: 

    

       # if(d.length != L):
        #    d.dis["P(x)"] = np.array( list(np.zeros(shape=int((L - d.length)/2))) + list(d.dis["P(x)"]) + list(np.zeros(shape=int((L - d.length)/2))) )
        #    d.nbin = data[0].nbin
         #   d.dis["ibin"] = data[0].dis["ibin"]

        binned, bin_mid, x_mid = ordered_binning(df=d,binsize=binsize)
        
        if type == 'linear':
            ax.set(xscale='linear',yscale='linear')
            ax.set(xlabel='$ibin$',ylabel='$P(x)$')
        if type == 'log':
            ax.set(xscale='linear',yscale='log')
            ax.set(xlabel='$x$',ylabel='$P(x)$')
        if type == 'gaus':
            ax.set(xscale='linear',yscale='log')
            ax.set(xlabel='$x^2$',ylabel='$P(x)$')
            bin_mid = bin_mid**2*np.sign(bin_mid)
        return binned,bin_mid
        ax.plot(bin_mid,binned)#,label=d.get_label(label))
    ax.legend()

    '''
    cols = len(data)
    title = ["nconf","length","gamma","nt","weight","t"]
    sup = "bin width=" + str(int(binsize)) + ", "
    for t in title:
        if t not in label: 
            sup += str(data[cols-1].get_label(t)) + ", "
    sup = sup[:-2]
    fig.suptitle(t=sup)
    '''



def ordered_binning(df, binsize, low_bound=None, high_bound=None):
    if low_bound == None: 
        low_bound = 0
        print("h1")
    if high_bound == None:
        high_bound = len(df.dis["P(|x|)"]) - 1
        print("h2")
    
    px = df.dis["P(|x|)"]  
    partitions = len(px[low_bound:high_bound + 1])/binsize
    # partitions = int(int(df.params["NBINS"])*2/binsize)
    print("num_partitions:", partitions)
    binned_dis = np.zeros(shape=int(partitions))
    midpoint_bins = np.zeros(shape=int(partitions))
    midpoint_x = np.zeros(shape=int(partitions))
    for n in range(int(partitions)):
        bin_start = low_bound + binsize*n
        binned_dis[n] = np.sum( px[bin_start:bin_start + binsize] )
        midpoint_bins[n] = df.dis["ibin"][bin_start + int(0.5*binsize)]
        midpoint_x[n] = df.dis["x"][bin_start + int(0.5*binsize)]

    return binned_dis, midpoint_bins, midpoint_x

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

        ax.plot(10**x_test, 10**series(x_test),marker='x',markersize=3,label=str(p.series))
        ax.legend()



def msd_fit_to_exp(d, interval,expr):
    x = d.avx["time"]
    y = d.avx["<r^2>"]

    nconfs = int(d.params["NSETS"]) * int(d.params["NWALKERS_PER_SET"])
    yerrors = np.sqrt( (d.avx["<r^4>"] - d.avx["<r^2>"] ** 2.0) / nconfs )

    bound_arr = x
    low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-7)[0][0] # find index where x==interval[0]
    high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-7)[0][0] + 1

    popt, pcov, infodict, errmsg, success = opt.curve_fit(expr, x[low_bound:high_bound], y[low_bound:high_bound], sigma=yerrors[low_bound:high_bound], full_output=True)
    # print(pcov)
    chi2 = sum ( infodict["fvec"] ** 2.0 )
    chi = np.sqrt ( chi2 / (len(y[low_bound:high_bound]) - len(popt)) )
    return yerrors, popt, pcov, chi2, chi



def plot_msd_with_fits(dfs, intvs, colors, markers, floating_exp_mode=True):
    fig,ax = plt.subplots(1,1)

    expr = lambda x, A, alpha: A * x ** alpha
    ax.set(xscale='log',yscale='log')
    markevery = 5
    ax.set(xlabel='t', ylabel='$ \langle r^2 \\rangle $')
    for i,df in enumerate(dfs):

        gamma = float(df.params["GAMMMA"])
        if (not floating_exp_mode): # one param mode 
            alpha = 2.0 - gamma 
            if(alpha < 4.0/3.0): 
                alpha = 4.0/3.0
            expr = lambda x, A: A * x ** alpha

        yerrors, popt, pcov, chi2, chi = msd_fit_to_exp(df, intvs[i], expr)
        label = f'$\\alpha={2.0 - gamma}$'
        t = df.avx["time"]
        ax.errorbar(t[::markevery], df.avx["<r^2>"][::markevery], yerrors[::markevery], label=label, color=colors[i], marker=markers[i], ls='none', markersize=6)
        ax.plot(t, expr(t, *popt), color=colors[i])
        print(f'gamma={gamma}, params={popt}, cov={np.sqrt(np.diag(pcov))}, chi2={chi2}, chi={chi}')
    ax.legend()


def log_fit2(x,y,interval,linear_error):
    bound_arr = x
    low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-7)[0][0] # find index where x==interval[0]
    high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-7)[0][0] + 1

    logx = np.log(x)
    logy = np.log(y)

    logerror = (linear_error / y)
    
    # weights w = 1/sigma, and here we need the logerror because we're fitting to the log data 
    series,stats = Poly.fit(logx[low_bound:high_bound], logy[low_bound:high_bound], deg=1, window=None, w=1/logerror[low_bound:high_bound], full=True)
    series = series.convert()

    chi2 = 0
    resid = sum(logy- series(logx))**2


    chi2 = sum( (logy[low_bound:high_bound] - series(logx[low_bound:high_bound]))**2 / logerror[low_bound:high_bound]**2  )
    dof = len(logy[low_bound:high_bound]) - 2 # 2 is number of fit params 
    chi = np.sqrt(chi2/dof)
    #chi2_lin = sum( (np.exp(logy[low_bound:high_bound]) - np.exp(series(logx[low_bound:high_bound])))**2 / (linear_error[low_bound:high_bound])**2  )

    lin_fit_over_intv = [ [ np.exp(x_i) for x_i in logx[low_bound:high_bound] ] , [ np.exp(series(x_i)) for x_i in logx[low_bound:high_bound] ] ]
    return lin_fit_over_intv, series, chi2, chi

def log_fit(ax,data_x,data_y,interval,plot_fit = True, weights=None):
    bound_arr = data_x
    low_bound = np.nonzero(np.fabs( bound_arr - interval[0]) < 1e-7)[0][0]
    high_bound = np.nonzero(np.fabs( bound_arr - interval[1]) < 1e-7)[0][0]

    x_test = np.linspace(np.log10(interval[0]),np.log10(interval[1]),100)
    # weights passed in as std. errors -> fit wants them as 1/std.err 
    series = Poly.fit(np.log10(data_x[low_bound:high_bound]), np.log10(data_y[low_bound:high_bound]), deg=1, window=None, w=1/weights[low_bound:high_bound]**2)
    series = series.convert()

    label = str(series) 
    if(len(weights)):
        
        # computes elements to be summed for xi^2 = sum ( (y_i - f(x_i))^2/sigma_i^2 ) 
        comps = (10**series(np.log10(data_x[low_bound:high_bound])) - 10**np.log10(data_y[low_bound:high_bound]))**2 / ((data_x[low_bound:high_bound]*weights[low_bound:high_bound])**2)
        print(comps)
        label += ", $\chi^2=$" + str( np.sum( comps )  ) 

    if (True):
        ax.plot(10**x_test, 10**series(x_test),markersize=4,marker='*',label=label)
        #ax.plot(x_test, series(x_test),marker='x',markersize=3)

    ax.legend()

    return series 

def gen_fit(ax,data_x,data_y,interval):
    bound_arr = data_x

    dif_low = abs(bound_arr[0] - interval[0])
    dif_high = abs(bound_arr[0] - interval[1])
    low_val = bound_arr[0]
    high_val = bound_arr[0]

    for x in data_x:
        if (abs(x - interval[0]) < dif_low):
            dif_low = abs(x - interval[0])
            low_val = x
        if (abs(x - interval[1]) < dif_high):
            dif_high = abs(x - interval[1])
            high_val = x

    print(low_val,high_val)
    low_bound = np.nonzero(np.fabs( bound_arr - low_val) < 1e-7)[0][0]
    high_bound = np.nonzero(np.fabs( bound_arr - high_val) < 1e-7)[0][0]

    x_test = np.linspace(interval[0],interval[1],100)
    series = Poly.fit(data_x[low_bound:high_bound], data_y[low_bound:high_bound], deg=1, window=None)
    series = series.convert()
    series = series.convert()

    ax.plot(x_test, series(x_test),markersize=3,label=series,linestyle='dashed')
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