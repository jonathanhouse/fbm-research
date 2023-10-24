import numpy as np
import random
import matplotlib.pyplot as plt

NT = 100
L = 10000000
nconf = 1000


sum_dis = np.zeros(L*2)



for c in range(nconf):
    dis = np.zeros(L*2)
    dis[L] = 1
    x = L
    for n in range(NT):
        r = random.random()

        if r > 0.5:
            x +=int(dis[x] - dis[x+1])
        if r < 0.5:
            x += int(dis[x-1] - dis[x])
        dis[x] += 1
    sum_dis[:] += dis[:]
    print(c)



plt.plot(sum_dis)
plt.show()