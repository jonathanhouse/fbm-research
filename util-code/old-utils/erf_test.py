import numpy as np
import matplotlib.pyplot as plt


def erfcc(x):
    z=abs(x)
    t=1.0/(1.0+0.5*z)
    erfcc_ret=t*np.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    if (x < 0.0): erfcc_ret=2.0-erfcc_ret
    return erfcc_ret


WALKS_PER_SET = 1
temp_xx = np.zeros(WALKS_PER_SET)
LEN_PER_BIN = 0.1
WALKER_WIDTH = 1.5
NBIN = 5000
L = 1000
conf_history = np.zeros(2*NBIN + 1)
WALKER_SIGMA = 3.0

for iwalker in range(WALKS_PER_SET):
    temp_xx[iwalker] = 499.9999019

for iwalker in range(WALKS_PER_SET):
    ibin = round(temp_xx[iwalker]/LEN_PER_BIN)
    print('ibin:', ibin)
    for j in range(-round(WALKER_WIDTH/LEN_PER_BIN), round(WALKER_WIDTH/LEN_PER_BIN) + 1, 1):
        print(ibin + j)
        if(abs(ibin + j) <= (NBIN)):
            #print(ibin + j)
            tmp_real = (erfcc(((ibin + j)*LEN_PER_BIN - temp_xx[iwalker])/(WALKER_SIGMA*np.sqrt(2.0))) - erfcc(((ibin + j + 1)*LEN_PER_BIN - temp_xx[iwalker])/(WALKER_SIGMA*np.sqrt(2.0))))/(2.0*LEN_PER_BIN)
            print(ibin + j, tmp_real)
            conf_history[ibin + j + 5000] = tmp_real + conf_history[ibin + j + 5000]


fig,ax = plt.subplots()
ax.plot(np.arange(-5000,5000+1,1), conf_history)
plt.show()
