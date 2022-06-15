import numpy as np
import random as r
import matplotlib.pyplot as plt

# Constants
mix     = 1.00              # Mix factor of new and old V
Rs      = 1.00              # Sample resistance
R0      = 1.00              # Load resistance
C0      = 1.20              # Coupling capacitance
Cp      = 0.50              # Sample capacitance
Vl_th   = 0.45              # Left sample V threshold
Vr_th   = 0.50              # Right sample V threshold
dV      = .005              # Prob. distr. width
dt      = .001*(1 + mix)    # Time-step
sim_len = 500               # Durantion of sim
iters   = 100               # Nuber of iterations

# Starting values
Vcl_0   = 0.0   # Initial left sample capacitor voltage
Vcr_0   = 0.0   # Initial right sample capacitor voltage

# Input
Vapp        = 1.0*np.ones(sim_len)      # Voltage input
time        = 1*np.arange(0, sim_len)   # Time-steps

# Counters
dtl = 0.0   # Time since last spike of left device
dtr = 0.0   # Time since last spike of right device

# Pre-factors
a   = (R0*(Rs*Cp + dt) - dt*Rs)/((Rs*Cp + dt)*R0)
b   = Rs*Cp/(Rs*Cp + dt)

# Lists
I   = []    # Total current
Il  = []    # Left branch current
Ir  = []    # Right branch current
I0  = []    # Coupling capacitor current
Vcl = []    # Left sample capacitor voltage
Vcr = []    # Right sample capacitor voltage
V0  = []    # Coupling capacitor voltage
tsl = []    # Spike timings of left device
tsr = []    # Spike timings of right device

# Compute current of left device (capacitor + sample)
def getIl(t):
    return (Vapp[t] + Vcl[-1])/R0    
 
# Compute current of right device (capacitor + sample) 
def getIr(t):
    return (Vapp[t] + Vcr[-1])/R0

# Compute voltage of coupling capacitor    
def getV0():
    return V0[-1] + I0[-1]*dt/C0

# Compute current of coupling capacitor   
def getI0():
    return - V0[-1]/(dt/C0 + 2*b*dt/(Cp*a)

# Compute voltage of left capacitor    
def getVcl(t, I0_new):
    return b/a*(Vcl[-1] + (Vapp[t]/R0 + I0[-1])*dt/Cp)

# Compute voltage of right capacitor    
def getVcr_new(t, I0_new):
    return b/a*(Vcr[-1] + (Vapp[t]/R0 - I0[-1])*dt/Cp)

# Compute probability of firing   
def getP(V, ts, V_th):
    if (V < V_th):
        return 1. - 1./(1. + ts/dt*np.exp((V - V_th)/dV))  
    else:
        e = np.exp(-(V - V_th)/dV) 
        return 1. - e/(e + ts/dt) 

# Solve equations of the circuit
def solveCircuit(t):
    I0.append(getI0())
    V0.append(getV0())
    Vcl.append(getVcl(t))
    Vcr.append(getVcr(t))
    Il.append(getIl(t))
    Ir.append(getIr(t))

# Mix the last and penultimate voltage values    
def mixSolutions(s = 'none'):
    if s == 'l':
        Vcl[-1] = (Vcl[-1] + mix*Vcl[-2])/(1 + mix)
    elif s == 'r':
        Vcr[-1] = (Vcr[-1] + mix*Vcr[-2])/(1 + mix)
    else:
        Vcl[-1] = (Vcl[-1] + mix*Vcl[-2])/(1 + mix)
        Vcr[-1] = (Vcr[-1] + mix*Vcr[-2])/(1 + mix)

for t in time:     
    if len(Vcl) == 0:
        Vcl.append(Vcl_0)
        Vcr.append(Vcr_0)
        V0.append(Vcl_0 - Vcr_0)
        
    solveCircuit()
    mixSolutions()
            
    rand = r.random()    
    if rand < getP(Vcl[-1], dtl, Vl_th): 
        tsl.append(t)
        dtl = 0.        
        Vcl[-1] = Vcl_0
    else:
        dtl = dtl + 1.
        
    if (rand < getP(Vcr[-1], dtr, Vr_th)):  
        tsr.append(t)
        dtr = 0.        
        Vcr[-1] = Vcr_0
    else:
        dtr = dtr + 1.
        
    t = t + 1

plt.figure()
fig, ax = plt.subplots(figsize=(15, 4))
ax.plot(time, Vcl[:sim_len], color='black', label='Vcl', marker='s')
ax.plot(time, Vcr[:sim_len], color='red', label='Vcr', marker='s')
plt.legend()
ax.set(xlabel='Time (arb. units)', ylabel='Voltage (arb. units)')
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
#ax.set_ylim(top=8)
#ax.yaxis.set_ticks(np.arange(0, 8, 1))
ax.tick_params(axis = 'both', which = 'both', labelsize = 16, length = 4)
fig.tight_layout()
plt.savefig('V(t).pdf')  

plt.figure()
fig, ax = plt.subplots(figsize=(15, 4))
ax.plot(time, Il, color='black', label='Il', marker='s')
ax.plot(time, Ir, color='red', label='Ir', marker='s')
plt.legend()
ax.set(xlabel='Time (arb. units)', ylabel='Output current (arb. units)')
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
#ax.set_ylim(top=8)
#ax.yaxis.set_ticks(np.arange(0, 8, 1))
ax.tick_params(axis = 'both', which = 'both', labelsize = 16, length = 4)
fig.tight_layout()
plt.savefig('I(t).pdf') 
    
    
    
    