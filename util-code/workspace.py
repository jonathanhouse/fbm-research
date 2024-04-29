from data_files import DataFile 
import matplotlib.pyplot as plt
import numpy as np
from plot import plot, msd_fit, plot_times, plot_binned, log_fit, gen_fit
from numpy.polynomial.polynomial import Polynomial as Poly
#import tikzplotlib as tz # i broke this kinda, need to use a different python intepretter in vscode than homebrew bc it doesn't have this 
# formerlly i had homebrew/bin as interpretter and used the python3 from the Frameworks/version/3.11 python
from data_files1 import DataFile1
import pickle

fig,ax = plt.subplots(1,1)

MID_BIN = 2*5000000

def suptitle_gen(x):
    title = "nbin/L=" + "5M" + "/" + "10M" + ", gamma=" + str(x.gamma) + ", nconf=" + str(x.nconf) + ", t=[0,2**26]" 
   
def suptitle_gen1(x):
    title = "nbin/L=" + "5M" + "/" + "10M" + ", gamma=" + str(x.params["GAMMMA"]) + ", nconf=" + str(x.params["NCONF"]) + ", NT=" + str(x.params["NT"]) 

    fig.suptitle(title)

run25 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
run26 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.99/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
#run94 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.9999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",'avx')
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
p1_test = DataFile("data/probabilistic force/paccept-weighted-walk/gamma=1.0/p_accept=0.9999/p-accept-conf/Foundry-32624.out",["ln(p_accept)="])


new_pa08 = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/nt=2**26/L=10M (nbin=5M)/nconf=20K",'avx')



run_g14 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=1.4/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g08 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.8/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",grab_dis=False)
run_g04 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.4/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g06 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128",'avx')
run_g1 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]",grab_dis=False)
run_g07 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.7/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2*26]",'avx')
run_g05 = DataFile("data/parallel walkers/procs as sets/linear force/gamma=0.5/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128",'avx')

#new_data = DataFile("data/probabilistic force/gamma=1.0/p_accept=0.8/test",'grad')

print(max(p1_test.markers),sum(p1_test.markers)/len(p1_test.markers),min(p1_test.markers))

run_fbm_94 = DataFile1("data/probabilistic force/gamma=1.0/p_accept=0.9999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",grab_dis=False)
run_fbm_95 = DataFile1("data/probabilistic force/gamma=1.0/p_accept=0.99999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",grab_dis=False)
run_fbm_93 = DataFile1("data/probabilistic force/gamma=1.0/p_accept=0.999/nt=2**26/L=10M (nbin=5M)/intv=[0,2**26]",grab_dis=False)


run_32_64 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=32*64",grab_dis=False)
run_16_128 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128",grab_dis=False)
run_64_32 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=64*32/asymmetric gradient (+1:32)",grab_dis=False)
run_128_16 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=128*16/asymmetric gradient/intv=[0,2**26]",grab_dis=False)
run_256_8 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=256*8/asymmetric symmetric/intv=[0,2**26]",grab_dis=False)
run_512_4 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=512*4/asymmetric gradient/intv=[0,2**26]",grab_dis=False)
run_1024_2 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=1024*2/asymmetric gradient/intv=[0,2**26]",grab_dis=False)
run_2048_1 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=2048*1/asymmetric gradient/intv=[0,2**26]",grab_dis=False)


avg_sum_test = DataFile1("data/probabilistic force/gamma=1.0/p_accept=0.9999/nt=2**26/L=10M (nbin=5M)/avg-sum-test",grab_dis=False)

asym_g1_fbmdet = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/fBM-grad-determiner/t_off=10/gamma=1.0",grab_dis=False)
asym_g4_fbmdet = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/fBM-grad-determiner/t_off=10/gamma=0.4",grab_dis=False)
asym_g7_fbmdet = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/fBM-grad-determiner/t_off=10/gamma=0.7/intv=[0,2**26]",grab_dis=False)
asym_g1_5050 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=1.0/intv=[0,2**26]",grab_dis=True)
asym_g4_5050 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=0.4",grab_dis=False)
asym_g7_5050 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=0.7",grab_dis=False)
symm_g4 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=symm/t_off=10/gamma=0.4",grab_dis=False)
symm_g7 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=symm/t_off=10/gamma=0.7",grab_dis=False)
symm_g1 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=symm/t_off=10/gamma=1.0/stepsig=1.0",grab_dis=False)

asym_g8_fbm = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/fBM-grad-determiner/t_off=10/gamma=0.8",grab_dis=False)
asym_g8_5050 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=0.8",grab_dis=False)
symm_g8 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=symm/t_off=10/gamma=0.8",grab_dis=False)

test_pa = DataFile1("data/probabilistic force/paccept-weighted-walk/gamma=1.0/p_accept=0.8",grab_dis=False)
test_pa1 = DataFile1("data/probabilistic force/paccept-weighted-walk/gamma=1.0/p_accept=0.9999",grab_dis=False)
test_pa2 = DataFile1("data/probabilistic force/paccept-weighted-walk/gamma=1.0/p_accept=0.99999",grab_dis=False)
test_p = DataFile1("data/probabilistic force/paccept-weighted-walk/gamma=1.0/p_accept=1.0/test",grab_dis=False)
test_indep =DataFile1("data/linear force/gamma=1.0/weight=0.0/nt=2**26/L=1.5M/nbin=750K",grab_dis=False)

test1 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]",grab_dis=False)
test05 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]/sigstep=0.50",grab_dis=False)
test025 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]/sigstep=0.25",grab_dis=False)
test75 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]/sigstep=0.75",grab_dis=False)
test01 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]/stepsig=0.1",grab_dis=False)
test001 = DataFile1("data/parallel walkers/procs as sets/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/intv=[0,2**26]/stepsig=0.01",grab_dis=False)
test0001 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=symm/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.0/stepsig=0.001",grab_dis=False)
test15 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=symm/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.0/stepsig=0.15",grab_dis=False)
test005 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=symm/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.0/stepsig=0.05",grab_dis=False)

log1 = DataFile("data/probabilistic force/paccept-weighted-walk/gamma=1.0/p_accept=0.9999/p-accept-conf-2/Foundry-65557.out",["ln(p_accept)="])

asym_g1_t1K = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=1.0/intv=[0,1K]",grab_dis=False)
asym_g1_t10K = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=1.0/intv=[0,10K]",grab_dis=False)
asym_g1_t1M = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=1.0/intv=[0,1M]",grab_dis=False)
asym_g1_tfull = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=1.0/intv=[0,2**26]",grab_dis=False)

gaus_1 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=gaus/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.0/intv=[0,2**11]",grab_dis=False)
gaus7 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=gaus/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=0.7/intv=[0,2**11]",grab_dis=False)
gaus13 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=gaus/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.3/intv=[0,2**11]",grab_dis=False)

gtest = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=unbinned/bin-form=gaus/t_off=10/gamma=1.0/intv=[0,2**11]",grab_dis=False)
gtest05 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=unbinned/bin-form=gaus/t_off=10/gamma=1.0/gaus-sig=0.5",grab_dis=False)
g04test = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=unbinned/bin-form=gaus/t_off=10/gamma=0.4/intv=[0,2**11]",grab_dis=False)


fs_test1 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.025/grad=unbinned/bin-form=gaus/t_off=10/gamma=1.0/intv=[0,2**12;fix1]",grab_dis=False)
fs_test2 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=gaus/weight=-0.025/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.0/intv=[0,2**12;fix1]",grab_dis=False)
fs_test3 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=symm/t_off=10/gamma=1.0/intv=[0,2**12]",grab_dis=False)
fs_test4 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=asym/50-50-determiner/t_off=10/gamma=1.0/intv=[0,2**26; force-steps]",grab_dis=False)

foff_test1 = DataFile1("data/parallel walkers/procs as sets/linear force/fbm_off/weight=-0.25/grad=unbinned/bin-form=gaus/t_off=10/gamma=1.0/intv=[0,2**12]",grab_dis=True)
fon_test2 = DataFile1("data/parallel walkers/procs as sets/linear force/grad=gaus/weight=-0.25/nt=2**26/L=10M (nbin=5M)/nconf=16*128/gamma=1.0/intv=[0,2**12]",grab_dis=True)

print("log1: ", max(log1.markers))

i = 0
for n in param_test.markers:
       if (n >= -3423070.9415688282):
              print(n)
              i += 1
print(i)


k = (2**26)*(32*64)*(1)*(1/64)*(1/32)


#y1= new_data.avx['<r>']
#y2 = new_data.avx['<r^2>']

#ax.errorbar(new_data.avx['t'],new_data.avx['<r>'],yerr=np.sqrt(y2-y1**2),label='<r^2>')
#ax.errorbar(new_data.avx['t'],new_data.avx['<r^2>'],yerr=np.sqrt(y2-y1**2))
#gen_fit(ax,np.sign(test_avx.dis['x'])*test_avx.dis['x']**2,np.log(test_avx.dis['P(x)']),interval=[1e10,0.5e11])
#ax.legend()
#suptitle_gen(new_data)

#plt.show()

#figx = pickle.load(open('msd_g8.fig.pickle', 'rb'))
#figx.show()

x_i = 't'
y_i = '<r>'
ax.set(xlabel='lnx',ylabel='ln(P(x))',xscale='linear',yscale='linear')
ax.set_title("16 * 128 walkers: grad_form=unbinned-gaus-sig=1.0; t_off=10")

print(test75.params)

tests = [foff_test1,fon_test2]
test_lab = ['fbm off bin_form=gaus,weight=0.25','fbm on bin_form=gaus,weight=0.25','asym-5050']

i = 0
for t in tests:
      
      ax.plot(t.full['ibin'],t.full['P(x)'],label=str(test_lab[i]))
       #log_fit(ax,t.full['time'],t.full['<r^2>'],interval=[1096,3687])
      i += 1
'''
log_fit(ax,tests[0].full["time"],tests[0].full["<r^2>"],interval=[1096,3687])
log_fit(ax,tests[1].full["time"],tests[1].full["<r^2>"],interval=[1096,3687])
log_fit(ax,tests[3].full["time"],tests[3].full["<r^2>"],interval=[1122548,60412582])
'''
suptitle_gen1(fs_test1)
ax.legend()
plt.show()


x = test75.dis["ibin"]
x2 = np.power(x,2)*np.sign(x)
ax.plot(x2,test75.dis["P(x)"],label='stepsig=0.75')

x = test05.dis["ibin"]
x2 = np.power(x,2)*np.sign(x)
ax.plot(x2,test05.dis['P(x)'],label='stepsig=0.5')

x = test025.dis["ibin"]
x2 = np.power(x,2)*np.sign(x)
ax.plot(x2,test025.dis['P(x)'],label='stepsig=0.25')

x = test01.dis["ibin"]
x2 = np.power(x,2)*np.sign(x)
ax.plot(x2,test15.dis["P(x)"],label='stepsig=0.15')

x = test15.dis["ibin"]
x2 = np.power(x,2)*np.sign(x)
ax.plot(x2,test01.dis["P(x)"],label='stepsig=0.1')

x = test005.dis["ibin"]
x2 = np.power(x,2)*np.sign(x)
ax.plot(x2,test005.dis["P(x)"],label='stepsig=0.05')
suptitle_gen1(test75)
ax.legend()
plt.show()


ax.plot(test75.avx["time"],test75.avx["<r^2>"],label='stepsig=0.75')
ax.plot(test05.avx["time"],test05.avx['<r^2>'],label='stepsig=0.5')
ax.plot(test025.avx["time"],test025.avx['<r^2>'],label='stepsig=0.25')
ax.plot(test15.avx['time'],test15.avx["<r^2>"],label='stepsig=0.15')
ax.plot(test01.avx["time"],test01.avx["<r^2>"],label='stepsig=0.1')
ax.plot(test005.avx["time"],test005.avx["<r^2>"],label='stepsig=0.05')
ax.plot(test001.avx['time'],test001.avx["<r^2>"],label='stepsig=0.01')
log_fit(ax,test001.avx['time'],test001.avx["<r^2>"],interval=[8980384,60412582])
#ax.plot(test0001.avx['time']/1.0e2,test0001.avx["<r^2>"]**(6.0e-1),label='stepsig=0.001')
ax.plot(test0001.avx['time'],test0001.avx["<r^2>"],label='stepsig=0.001')
log_fit(ax,test0001.avx['time'],test0001.avx["<r^2>"],interval=[8980384,60412582])
log_fit(ax,test0001.avx['time'],test0001.avx["<r^2>"],interval=[1,1096])

suptitle_gen1(test01)
ax.legend()
plt.show()




#log_fit(ax,gaus_1.avx["time"],gaus_1.avx["<r^2>"],interval=[114,1843])
ax.plot(gaus_1.dis["ibin"],gaus_1.dis["P(x)"],label='gamma=1.0')

w = gaus7.dis
x = w['ibin']
y = w['P(x)']
#log_fit(ax,x,y,interval=[114,1843])
ax.plot(x,y,label='gamma=0.7')

w = gaus13.dis
x = w['ibin']
y = w['P(x)']
#log_fit(ax,x,y,interval=[1096,1843])
ax.plot(x,y,label='gamma=1.3')

ax.legend()
suptitle_gen1(gaus_1)
plt.show()



w = asym_g1_tfull.dis
x = w["ibin"]
y2 = w["P(x)"]
ax.errorbar(np.sign(x)*np.power(x,2)*(1.0e4/2**26)**(4./3),y2*(1.0e4/2**26)**(1./3)*(1.0e3)*6.7,label='asym50505;t=[0,2**26]' )

w = asym_g1_t1M.dis
x = w["ibin"]
y2 = w["P(x)"]
ax.errorbar(np.sign(x)*np.power(x,2)*(1.0e4/1.0e6)**(4./3),y2*(1.0e4/1.0e6)**(1./3)*(1.0e1)*10,label='asym50505;t=[0,1M]' )

w = asym_g1_t10K.dis
x = w["ibin"]
y2 = w["P(x)"]
ax.errorbar(np.sign(x)*np.power(x,2),y2,label='asym50505;t=[0,10K]' )

suptitle_gen1(asym_g1_tfull)
ax.legend()
plt.show()
'''
w = asym_g1_t10K.dis
x = w["ibin"]
y2 = w["P(x)"]
ax.errorbar(np.sign(x)*np.power(x,2)/10,(10)**(4/3.)*y2,label='asym5050;t=[0,10K]' )

ax.set_xlim(-10000,10000)
'''

'''
w = asym_g1_t1M.dis
x = w["ibin"]
y2 = w["P(x)"]
ax.errorbar(np.sign(x)*np.power(x,2),y2,label='asym5050;t=[0,1M]' )
'''

#w = asym_g1_5050.dis
#x = w["ibin"]
#y2 = w["P(x)"]
#ax.errorbar(np.sign(x)*np.power(x,2),y2,label='symm;t=[0,2**26]')


'''
'''


ax.legend()

suptitle_gen1(asym_g1_5050)
#pickle.dump(fig,open('msd_g8.fig.pickle', 'wb'))

plt.show()


'''
w = run_fbm_off.dis
ax.plot(np.sign(w['ibin'])*np.log(w['ibin']),np.log(w['P(x)']),label='gamma=1.0')


w = run_fbm_off7.dis
ax.plot(np.sign(w['ibin'])*np.log(w['ibin']),np.log(w['P(x)']),label='gamma=0.7')

w = run_fbm_off4.dis
ax.plot(np.sign(w['ibin'])*np.log(w['ibin']),np.log(w['P(x)']),label='gamma=0.4')

ax.legend()



suptitle_gen(run_fbm_off)

plt.show()

'''
'''
w = run_g08.dis
x = w['ibin']
y = w['P(x)']
ax.plot(x,y,label='gamma=0.8')

w = run_g1.dis
x = w['ibin']
y = w['P(x)']
ax.plot(x,y,label='gamma=1.0')
'''




