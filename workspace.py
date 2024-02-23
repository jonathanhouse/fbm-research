from data_files import DataFile 
import matplotlib.pyplot as plt 
import numpy as np
from plot import plot, msd_fit, plot_times, plot_binned, log_fit, gen_fit
from numpy.polynomial.polynomial import Polynomial as Poly
import tikzplotlib as tz 


fig,ax = plt.subplots(1,1)

def suptitle_gen(x):
    title = "nbin/L=" + "5M" + "/" + "10M" + ", gamma=" + str(x.gamma) + ", nconf=" + str(x.nconf) + ", t=[0,2**26]" 
   

    fig.suptitle(title)

run25 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run26 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.99/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run94 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.9999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run01 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.1/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run05 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.5/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')

run_wavx = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/nt=2**26/L=10M (nbin=5M)/weighted avg/intv=[0,2**26]",'avx')

dist1 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=64*32/asymmetric gradient (+1:32)",'avx')
run1walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=2048*1/asymmetric gradient/intv=[0,2**26]",'avx')
run2walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=1024*2/asymmetric gradient/intv=[0,2**26]",'avx')
run4walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=512*4/asymmetric gradient/intv=[0,2**26]",'avx')
run8walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=256*8/asymmetric symmetric/intv=[0,2**26]",'avx')
run16walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=128*16/asymmetric gradient/intv=[0,2**26]",'avx')
run64walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=32*64",'avx')
run128walker = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128",'avx')
run1 = DataFile("data/linear force/gamma=1.0/weight=0.0/nt=2**26/L=1.5M/nbin=750K",'avx')

run2walker_rng = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=1024*2/asymmetric gradient/intv=[0,2**26]/IRINT=2048",'avx')

run59 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.99999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run69 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.999999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run79 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.9999999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run59_norm = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.99999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]/normalized",'avx')
run69_norm=  DataFile("data/probabilistic force/gamma=1.0/p_accept=0.999999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]/normalized",'avx')
run79_norm = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.9999999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]/normalized",'avx')
run1_norm = DataFile("data/probabilistic force/gamma=1.0/p_accept=1.0/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]/normalized",'avx')

wavg8 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/nt=2**26/L=10M (nbin=5M)/weighted avg (w norm)",'avx')
wavg92 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.99/nt=2**26/L=10M (nbin=5M)/weighted avg (w norm)",'avx')
wavg94 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.9999/nt=2**26/L=10M (nbin=5M)/weighted avg (w norm)",'avx')
wavg5 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.5/nt=2**26/L=10M (nbin=5M)/weighted avg (w norm)",'avx')

avg94 = DataFile("figures/winter break/final p_accept/p_accept=0.99/Foundry-2861961.out",["log(p_accept_conf)="])
avg9 = DataFile("figures/winter break/final p_accept/p_accept=0.9/Foundry-2861965.out",["log(p_accept_conf)="])
avg8 = DataFile("figures/winter break/final p_accept/p_accept=0.8/Foundry-2861966.out",["log(p_accept_conf)="])
avg5 = DataFile("figures/winter break/final p_accept/p_accept=0.5/Foundry-2861967.out",["log(p_accept_conf)="])
avg01 = DataFile("figures/winter break/final p_accept/p_accept=0.1/Foundry-2861971.out",["log(p_accept_conf)="])
avg03 = DataFile("figures/winter break/final p_accept/p_accept=0.3/Foundry-2861972.out",["log(p_accept_conf)="])
avg001 = DataFile("figures/winter break/final p_accept/p_accept=0.01/Foundry-2861973.out",["log(p_accept_conf)="])
avg0001 = DataFile("figures/winter break/final p_accept/p_accept=0.001/Foundry-2861974.out",["log(p_accept_conf)="])


avg6_01 = DataFile("figures/winter break/final p_accept/gamma=0.6/p_accept=0.1/Foundry-2862070.out",["log(p_accept)="])
avg6_05 = DataFile("figures/winter break/final p_accept/gamma=0.6/p_accept=0.5/Foundry-2862059.out",["log(p_accept)="])
avg6_08 = DataFile("figures/winter break/final p_accept/gamma=0.6/p_accept=0.8/Foundry-2862060.out",["log(p_accept_conf)="])
avg6_09 = DataFile("figures/winter break/final p_accept/gamma=0.6/p_accept=0.9/Foundry-2862065.out",["log(p_accept_conf)="])

param_test = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/p_accept_conf/Foundry-2869514.out",["log(p_accept)="])

new_pa08 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/nt=2**26/L=10M (nbin=5M)/nconf=20K",'avx')



run_g14 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=1.4/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g08 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.8/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g04 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.4/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g06 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128",'avx')
run_g1 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]",'avx')
run_g07 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.7/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g05 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.5/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128",'avx')

run_fbm_off = DataFile("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/t_off=10",'avx')
run_fbm_off4 = DataFile("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/t_off=10 (gamma=0.4)",'avx')


test_avx = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/test",'avx')

print(max(param_test.markers),sum(param_test.markers)/len(param_test.markers),min(param_test.markers))


i = 0
for n in param_test.markers:
       if (n >= -3423070.9415688282):
              print(n)
              i += 1
print(i)

x_i = 't'
y_i = '<r^2>'
ax.set(xlabel=x_i,ylabel=y_i)

ax.set(xscale='log',yscale='log')
ax.set(title='(16 sets, 128 walkers/set): mean-squared displacement with error-bars')
k = (2**26)*(32*64)*(1)*(1/64)*(1/32)

#ax.plot(x.cor['t'],(x.cor['<r^2>']))##

#                                                       ax.plot(test_avx.avx['t'],test_avx.avx['<r^2>'])

#plt.show()

dis25 = run25.avx
dis26 = run26.avx
data94 = run94.avx
data01 = run01.avx

data05 = run05.avx
#data_dist = dist1.dis
data1 = run1.avx

fit = np.power(10,-0.29153843)*np.power(dis25['t'],1.3354501)

#print((len(fit)))
tests = [data01,data05,dis25,dis26,data94]

i = 0
a = np.zeros(len(tests))
for test in tests:
       for t in range(0,99):
              x2_test = test['<r^2>'][t]
              x2_bm = data1['<r^2>'][t]

              cond = abs(x2_test-x2_bm)/(x2_test+x2_bm) > 0.2
              higher = (np.log(fit[t]) > np.log(x2_bm)) and (np.log(x2_test) > np.log(x2_bm)) and (np.log(x2_test) < np.log(fit[t]))

              if(cond):
                     #print(test['t'][t])
                     a[i] = test['t'][t]
                     break 

              if(t==98):
                     a[i] = 0
                     print("no pass")
       i += 1

#print(a)



p_accept = np.array([0.1,0.5,0.8,0.99,0.9999])

p_accept = np.exp(np.exp(-1/(1-p_accept)))



y_val = np.linspace(10e1,10e9,100)

for i in a:
       t_range = i*np.ones(100)
       y_val = np.linspace(1,10e9,100)
       #ax.plot(t_range,y_val,label='t='+str(i))

'''
ax.plot(t_range,y_val)
ax.plot(data01[x_i],data01[y_i], label='p_accept=0.1')
ax.plot(data05[x_i],data05[y_i], label='p_accept=0.5')
ax.plot(dis25[x_i],dis25[y_i], label='p_accept=0.8')
##ax.plot(dis25[x_i],fit, label='fit')
ax.plot(dis26[x_i],(dis26[y_i]), label='p_accept=0.99')
ax.plot(data94[x_i],(data94[y_i]), label='p_accept=0.9999')
ax.plot(data1[x_i],data1[y_i], label='BM')
'''
#ax.plot(data_dist[x_i],k*data_dist[y_i])
#ax.plot(data01[x_i],abs(fit-data01[y_i]), label='p_accept=0.1')
#ax.plot(data05[x_i],abs(fit-data05[y_i]), label='p_accept=0.5')

'''
                                         
ax.plot(run25.avx[x_i],run25.avx[y_i],label='unweighted avg')
ax.plot(run_wavx.avx[x_i],run_wavx.avx[y_i],label='weighted avg')
'''



'''
avg1 = sum(avg94.markers)/len(avg94.markers)
avg2 = sum(avg9.markers)/len(avg9.markers)
avg3  = sum(avg8.markers)/len(avg8.markers)
avg4 = sum(avg5.markers)/len(avg5.markers)
avg05 = sum(avg01.markers)/len(avg01.markers)
avg6 = sum(avg03.markers)/len(avg03.markers)
avg7 = sum(avg001.markers)/len(avg001.markers)
avg8a = sum(avg0001.markers)/len(avg0001.markers)

print(avg6_08.markers)

avg1_8 = sum(avg6_08.markers)/len(avg6_08.markers)
avg1_9 =  sum(avg6_09.markers)/len(avg6_09.markers)
avg1_1 =  sum(avg6_01.markers)/len(avg6_01.markers)
avg1_5 =  sum(avg6_05.markers)/len(avg6_05.markers)


a = [avg8a, avg7, avg05, avg6, avg4, avg3, avg2, avg1]
p = [0.001, 0.01, 0.1 ,0.3 , 0.5,0.8,0.9,0.99]

a1 = [avg1_1,avg1_5,avg1_8,avg1_9]
p1 = [0.1,0.5,0.8,0.9]

ax.scatter(np.log(p),a,label='gamma=1.0')
gen_fit(ax=ax,data_x=np.log(p),data_y=a,interval=[np.log(0.001),np.log(0.99)])

ax.scatter(np.log(p1),a1,label='gamma=0.6')

c1 = np.polyfit(np.log(p1),a1,1)
s1 = np.poly1d(c1)

x1 = np.linspace(np.log(min(p)),np.log(max(p)))
'''




#ax.plot(x1,s1(x1),label=s1)



'''
ax.plot(run05.avx[x_i],run05.avx[y_i],label='unnorm: paccept=0.5')
ax.plot(run25.avx[x_i],run25.avx[y_i],label='unnorm: paccept=0.8')
ax.plot(run26.avx[x_i],run26.avx[y_i],label='unnorm: paccept=0.99')
ax.plot(run94.avx[x_i],run94.avx[y_i],label='unnorm: paccept=0.99')

#ax.plot(run59_norm.avx[x_i],run59_norm.avx[y_i],label='p_accept=0.99999')
ax.plot(wavg5.avx[x_i],wavg5.avx[y_i],label='p_accept=0.5')
ax.plot(wavg8.avx[x_i],wavg8.avx[y_i],label='p_accept=0.8')
ax.plot(wavg92.avx[x_i],wavg92.avx[y_i],label='p_accept=0.99')
ax.plot(wavg94.avx[x_i],wavg94.avx[y_i],label='p_accept=0.9999')
'''


#msd_fit(ax=ax,data=[run1_norm],interval=[35921537,60412582])

'''
ax.plot(run1walker.avx[x_i],run1walker.avx[y_i],label='2048*1')
ax.plot(run2walker.avx[x_i],run2walker.avx[y_i],label='1024*2')
ax.plot(run4walker.avx[x_i],run4walker.avx[y_i],label='512*4')
ax.plot(run8walker.avx[x_i],run8walker.avx[y_i],label='256*8')
ax.plot(run16walker.avx[x_i],run16walker.avx[y_i],label='128*16')
ax.plot(dist1.avx[x_i],dist1.avx[y_i],label='64*32')
ax.plot(run64walker.avx[x_i],run64walker.avx[y_i],label='32*64')
'''


w = run_g04.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=0.4; ' + str(s))

w = run_g05.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=0.5; ' + str(s))

w = run_g06.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=0.6; ' + str(s))


w = run_g07.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=0.7; ' + str(s))

w = run_g08.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=0.8; ' + str(s))

w = run_g1.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=1.0; ' + str(s))

w = run_g14.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[1096,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='gamma=1.4; ' + str(s))

w = run_fbm_off.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[10429,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='force-only (gamma_init=1.0,t_off=10); ' + str(s))

w = run_fbm_off4.avx
x = w['t']
y1 = w['<r>']
y2 = w['<r^2>']
s = log_fit(ax,x,y2,interval=[10429,60412582])
ax.errorbar(x,y2,yerr=np.sqrt(y2-y1**2),label='force-only (gamma_init=0.4,t_off=10); ' + str(s))

ax.legend()

suptitle_gen(run_g06)

plt.show()



#plt.show()


w = run1walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]


ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='2048*1')

w = run2walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='1024*2')

w = run4walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='512*4')

w = run8walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='256*8')

w = run16walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='128*16')

w = dist1.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='64*32')

w = run64walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='32*64')

w = run128walker.avx
x2 = w["<r^2>"]
x1 = w["<r>"]
ax.errorbar(w[x_i],w[y_i],yerr=np.sqrt(x2-x1**2),label='16*128')

suptitle_gen(run1walker)



#ax.plot(new_pa08.avx[x_i],new_pa08.avx[y_i])

ax.legend()

#ax.ploax.plot(pos['t'],pos['<r^2>'])t(dis['ibin'],dis['P(x)'])
#ax.plot(pos['t'],pos['<r^2>'],label='nconf=64*128')
#ax.plot(pos1['t'],pos1['<r^2>'],label='nconf=64*64')

#ax.plot(pos3['t'],pos3['<r^2>'],label='nconf=10*64,non-parallel')

#log_fit(ax,dis25['t'],dis25['<r^2>'],interval=[35921537,60412582])


#ax.legend()
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
