import numpy as np
import random as r
import matplotlib.pyplot as plt

# Constants
Rs      = 0.3   # Sample resistance
R0      = 0.1   # Load resistance
C0      = 0.0   # Coupling capacitance
Cp      = 1.5   # Sample capacitance
Vl_th   = 10.0  # Left sample V threshold
Vr_th   = 10.0  # Right sample V threshold
dV      = 0.1   # Prob. distr. width
dt      = 1.0   # Time-step
sim_len = 100   # Durantion of sim

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
a   = R0*(Rs*Cp + dt)/(R0*Rs*Cp + dt*(R0 - Rs))
b   = Rs*Cp/(2*R0*(Rs*Cp + dt))
c   = Rs*Cp/(-Rs*Cp + dt)

# Lists
Il  = []    # Left branch current
Ir  = []    # Right branch current
I0  = []    # Coupling capacitor current
Vcl = []    # Left sample capacitor voltage
Vcr = []    # Right sample capacitor voltage
V0  = []    # Coupling capacitor voltage

def getIl_new(I_new, Vcl_old, Vcr_old, V0_old):
    return a*(I_new/2. + b*(-I_new*dt/Cp + Vcl_old - Vcr_old - 2*V0_old))
    
def getIr_new(I_new, Vcl_old, Vcr_old, V0_old):
    return a*(I_new/2. - b*(I_new*dt/Cp + Vcl_old - Vcr_old - 2*V0_old))
    
def getV0_new(I_new, Il_new, Vcl_old, Vcr_old, V0_old): 
    return -c*(dt/Cp*(2.*Il_new - I_new) + Vcl_old - Vcr_old - 2*V0_old)
    
def getI0_new(V0_new, V0_old):
    return (V0_new - V0_old)*C0
    
def getVcl_new(Il_new, I0_new, Vcl_old):
    return Vcl_0 if Cp == 0 else Rs*Cp/(Rs*Cp + dt)*(Vcl_old + (Il_new + I0_new)/Cp)
    
def getVcr_new(Ir_new, I0_new, Vcr_old):
    return Vcr_0 if Cp == 0 else Rs*Cp/(Rs*Cp + dt)*(Vcr_old + (Ir_new - I0_new)/Cp)
    
def getP(V, ts, V_th, dV):
    if (V - V_th) < 0.:
        return 1. - 1./(1. + ts*np.exp((V - V_th)/dV))
    else:
        e = np.exp(-(V - V_th)/dV)
        return 1. - e/(e + ts)  
    
for I_new in I:
    if len(Vcl) == 0:
        Vcl.append(Vcl_0)
        Vcr.append(Vcr_0)
        V0.append(V0_0)
        
    Il.append(getIl_new(I_new, Vcl[-1], Vcr[-1], V0[-1]))
    Ir.append(getIr_new(I_new, Vcl[-1], Vcr[-1], V0[-1]))
    V0.append(getV0_new(I_new, Il[-1], Vcl[-1], Vcr[-1], V0[-1]))
    I0.append(getI0_new(V0[-1], V0[-2]))
    Vcl.append(getVcl_new(Il[-1], I0[-1], Vcl[-1]))
    Vcr.append(getVcr_new(Ir[-1], I0[-1], Vcr[-1]))
    
    rand = r.random()
    
    if (rand < getP(Vcl[-1], tls, Vl_th, dV)):
        Vcl[-1] = Vcl_0
        tls = 0.
    else:
        tls = tls + 1.
        
    if (rand < getP(Vcr[-1], trs, Vr_th, dV)):
        Vcr[-1] = Vcr_0
        trs = 0.
    else:
        trs = trs + 1.

plt.figure()
fig, ax = plt.subplots(figsize=(6, 5))
ax.plot(t, Vcl[:sim_len], color='blue', marker='.', label='Vcl')
ax.plot(t, Vcr[:sim_len], color='red', marker='.', label='Vcr')
ax.set(xlabel='Time (arb. units)', ylabel='Voltage (arb. units)')
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
#ax.set_ylim(top=8)
#ax.yaxis.set_ticks(np.arange(0, 8, 1))
ax.tick_params(axis = 'both', which = 'both', labelsize = 16, length = 4)
fig.tight_layout()
plt.savefig('Vcl(t).pdf')  
    
    
    
    