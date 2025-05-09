import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


NT = 2**26
L = 1500000
NBIN = 750000


NTSTART =2**26-1000
NTEND = 2**26
NWATCH = NTEND - NTSTART + 1

IT = 0
XXNF = IT + 2
XXNI = XXNF + 2
FBMSTEP = XXNI + 2
FORCESTEP = FBMSTEP + 2
NEIGHBORDIS = FORCESTEP + 1

pause = False

#path = '../data/linear force/gamma=0.6/weight=-0.25/nt=2**26/L=10M/gradient_dx/dx=10/data.out'
path = '../data/linear force/gamma=1.0/weight=-0.25/nt=2**26/L=1.5M/new output/mixed gradient /nbin=750K/Foundry-2772471.out'
pos = np.zeros(NT)
dis = np.zeros(NBIN*2 + 1)
lby2 = int(L/2)

file_read = open(path,'r').readlines()
N = len(file_read)
print("processed file of size",N)

for i in range(N):
    line = file_read[i].split()
    if len(line) > 0:
        if line[0] == 'dis.':
            offset = i + 1
            break
#from data, get list of walker's x pos's
#

fig, ax = plt.subplots(1,1)
walker, = ax.plot([0],[0],'ro',markersize=5)
dis_plot, = ax.plot([],[])
grad_plot, = ax.plot([],[])

error = False
for t in range(0,NTSTART):

    pos[t] = file_read[t+offset].split()[XXNI] # t-1 is needed to grab the first line to fill t=1 spot of pos 


    #in_range = (pos[t] < pos[NTSTART] + 1000*NWATCH) and (pos[t] > pos[NTSTART] - 1000*NWATCH)
    #if(t < NTSTART and in_range): 

    x = pos[t]
    ibin = round(x*(NBIN/lby2))
    dis[ibin + NBIN] += 1
 # first error at t=2947
    neighbors = np.array(file_read[t+offset].split()[NEIGHBORDIS][1:-1].split(','),dtype=int)

    if (abs(ibin) != NBIN):
        if( dis[ibin-1 + NBIN] != neighbors[0] or 
            dis[ibin + NBIN] != neighbors[1] or 
            dis[ibin+1 + NBIN] != neighbors[2]):
                if error == False:
                    print("error at bin ", ibin, " at time ", t)
                    error = True

    else:
        if ( dis[ibin + NBIN] != neighbors[1] ):
            if error == False:
                print("error at bin ", ibin, " at time ", t, " at wall ", np.sign(ibin)*NBIN)
                error = True



print("loaded data into position vector")
ax.set(xlim=[-lby2,lby2],ylim=[0,2*max(dis)])

W = 5
weight = -0.25
def animate(t):
        
    t_curr = t + NTSTART
    t_curr_data = file_read[t_curr+offset].split()

    pos[t_curr] = t_curr_data[XXNI]
    x_curr =  pos[t_curr]

    ibin = round(x_curr*(NBIN/lby2))
    dis[ibin + NBIN] += 1
    
    fbm_step = float(t_curr_data[FBMSTEP])
    force_step = float(t_curr_data[FBMSTEP + 2])
    ax.set(title='t=' + str(t+NTSTART) +" (" + str(t_curr_data[IT]) + "\n" + str(t_curr_data[XXNF]) + " = " + str(x_curr) + " + " + str(fbm_step) + " + " + str(force_step))
    walker.set(color='red')
    if isinstance(fbm_step, str):     
        print(fbm_step)

    if fbm_step < 0:
        walker.set(color='blue')

    #y_mean = 0
    #for j in range(x - W, x + W + 1):  
    #    y_mean += dis[j + lby2]
    #intc = y_mean/(2*W + 1.) - force_step/weight * x
    #window = np.arange(x - W, x + W)


    walker.set_data(ibin*(lby2/NBIN),dis[ibin+NBIN])
    #grad_plot.set_data(window, window*force_step/weight + intc)
    dis_plot.set_data(np.linspace(-lby2,lby2,2*NBIN+1),dis)

    #

def on_click(event):
    global pause
    pause ^= True


'''
t_step = NTSTART
def next_frame():
    t_step += 1
    if t_step == NTEND:
        t_step = NTSTART

#keyboard.on_press_key("right arrow", lambda _:next_frame())
'''

fig.canvas.mpl_connect('button_press_event',on_click)
anim = FuncAnimation(fig,animate,frames=NWATCH,interval=1000)



plt.show()
