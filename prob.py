import numpy as np
import matplotlib.pyplot as plt

#Constants
V_th    = 0.50          # V threshold
dV      = .005          # Prob. distr. width
dt      = .001          # Time-step for traces
sim_len = 0.10     # Durantion of sim 

# Input
Vs      = np.arange(V_th-3*dV, V_th, dV) # Sample voltage
time    = np.arange(0, sim_len, dt)     # Time-steps

# Lists
P = []      # Probability
cum = []    # Cumulative

def getP(V, ts):
    return 1. - 1./(1. + ts/dt*np.exp((V - V_th)/dV))  

fig, ax = plt.subplots(figsize=(8, 4))

for v in Vs:
    cum = []
    P = getP(v, time)
    c = 0
    for p in P:
        cum.append(c)
        c += p
    ax.plot(time, P, label='V='+str(np.round(v,3))+' (arb. units)')
    
ax.set(xlabel='Time (arb. units)', ylabel='Probability')
size = 24 
ax.xaxis.label.set_size(size)
ax.yaxis.label.set_size(size)
ax.set_ylim(bottom=0.0)
ax.set_ylim(top=1.0)
ax.set_xlim(left=-0.0)
ax.set_xlim(right=sim_len)
plt.legend( prop={'size': 16})
ax.tick_params(axis = 'both', which = 'both', labelsize = size, length = int(size/6))
fig.tight_layout()
plt.savefig('./figures/prob.pdf') 
plt.close()