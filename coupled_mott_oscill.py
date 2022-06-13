import numpy as np
import random as r
import matplotlib.pyplot as plt

# Constants
Rs      = 1.00  # Sample resistance
R0      = 1.00  # Load resistance
C0      = 0.50  # Coupling capacitance
Cp      = 0.50  # Sample capacitance
Vl_th   = 0.50  # Left sample V threshold
Vr_th   = 0.45  # Right sample V threshold
dV      = 0.01  # Prob. distr. width
dt      = 0.0001  # Time-step
sim_len = 20000  # Durantion of sim

# Starting values
Vcl_0   = 0.0   # Initial left sample capacitor voltage
Vcr_0   = 0.0   # Initial right sample capacitor voltage
V0_0    = 0.0   # Initial coupling capacitor voltage

# Input
I   = 1.0*np.ones(sim_len)      # Input current
t   = 1.0*np.arange(0, sim_len) # Time-steps

# Counters
tls = 0.0   # Time since last spike of left device
trs = 0.0   # Time since last spike of right device

# Pre-factors
a   = R0*(3*Rs*Cp + dt)/(3*R0*Rs*Cp + dt*(R0 - Rs))
b   = Rs*Cp/(2*R0*(3*Rs*Cp + dt))
c   = Rs*Cp/(3*Rs*Cp + dt)

# Lists
Il  = []    # Left branch current
Ir  = []    # Right branch current
I0  = []    # Coupling capacitor current
Vcl = []    # Left sample capacitor voltage
Vcr = []    # Right sample capacitor voltage
V0  = []    # Coupling capacitor voltage

def getIl_new(I_new, Vcl_old, Vcr_old, V0_old):
    return a*(I_new/2. - b*(I_new*dt/Cp - Vcl_old + Vcr_old + 2*V0_old))
    
def getIr_new(I_new, Vcl_old, Vcr_old, V0_old):
    return a*(I_new/2. + b*(-I_new*dt/Cp - Vcl_old + Vcr_old + 2*V0_old))
    
def getV0_new(I_new, Il_new, Ir_new, Vcl_old, Vcr_old, V0_old): 
    return 0 if C0 == 0 else c*(dt/Cp*(Ir_new - Il_new) + Vcr_old - Vcl_old + 2*V0_old)
    
def getI0_new(V0_new, V0_old):
    return (V0_new - V0_old)*C0/dt
    
def getVcl_new(Il_new, I0_new, Vcl_old):
    return Vcl_0 if Cp == 0 else Rs*Cp/(Rs*Cp + dt)*(Vcl_old + (Il_new + I0_new)*dt/Cp)
    
def getVcr_new(Ir_new, I0_new, Vcr_old):
    return Vcr_0 if Cp == 0 else Rs*Cp/(Rs*Cp + dt)*(Vcr_old + (Ir_new - I0_new)*dt/Cp)
    
def getP(V, ts, V_th, dV):
    if (V < V_th):
        return 1. - 1./(1. + ts/dt*np.exp((V - V_th)/dV))  
    else:
        e = np.exp(-(V - V_th)/dV) 
        return 1. - e/(e + ts/dt) 

def solveCircuit(Vcl_old, Vcr_old):
    Il.append(getIl_new(I_new, Vcl_old, Vcr_old, V0[-1]))
    Ir.append(getIr_new(I_new, Vcl_old, Vcr_old, V0[-1]))
    V0.append(getV0_new(I_new, Il[-1], Ir[-1], Vcl_old, Vcr_old, V0[-1]))
    I0.append(getI0_new(V0[-1], V0[-2]))
    Vcl.append(getVcl_new(Il[-1], I0[-1], Vcl_old))
    Vcr.append(getVcr_new(Ir[-1], I0[-1], Vcr_old))

for I_new in I:
    if len(Vcl) == 0:
        Vcl.append(Vcl_0)
        Vcr.append(Vcr_0)
        V0.append(V0_0)
    
    if tls == 0.:
        solveCircuit(Vcl_0, Vcr[-1])
    elif trs == 0.:
        solveCircuit(Vcl[-1], Vcr_0)
    else:
        solveCircuit(Vcl[-1], Vcr[-1])
            
    rand = r.random()    
    if rand < getP(Vcl[-1], tls, Vl_th, dV):        
        tls = 0. 
        print('left')
        continue
    elif tls == 0.:
        Vcl[-1] = Vcl_0
        tls = tls + 1.
    else:
        tls = tls + 1.
        
    if (rand < getP(Vcr[-1], trs, Vr_th, dV)):         
        trs = 0.
        print('right')
    elif trs == 0.:
        Vcr[-1] = Vcr_0
        trs = trs + 1.
    else:
        trs = trs + 1.

plt.figure()
fig, ax = plt.subplots(figsize=(15, 4))
ax.plot(t, Vcl[:sim_len], color='blue', label='Vcl', marker='.')
ax.plot(t, Vcr[:sim_len], color='red', label='Vcr', marker='.')
plt.legend()
ax.set(xlabel='Time (arb. units)', ylabel='Voltage (arb. units)')
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
#ax.set_ylim(top=8)
#ax.yaxis.set_ticks(np.arange(0, 8, 1))
ax.tick_params(axis = 'both', which = 'both', labelsize = 16, length = 4)
fig.tight_layout()
plt.savefig('Vcl(t).pdf')  
    
    
    
    