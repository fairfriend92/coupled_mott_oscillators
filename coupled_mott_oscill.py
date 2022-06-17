import numpy as np
import random as r
import matplotlib.pyplot as plt

# Constants
Rs      = 1.00          # Sample resistance
R0      = 1.00          # Load resistance
Cp      = 0.50          # Sample capacitance
Vl_th   = 0.45          # Left sample V threshold
Vr_th   = 0.50          # Right sample V threshold
dV      = .005          # Prob. distr. width
dt      = .01           # Time-step
sim_len = int(100/dt)   # Durantion of sim
tRef    = 5             # Refractory period of device

# Parameters
C0Array  = [0., 0.5, 10.]    # Coupling capacitance

# Starting values
Vcl_0   = 0.0   # Initial left sample capacitor voltage
Vcr_0   = 0.0   # Initial right sample capacitor voltage

# Input
Vapp    = 1.0*np.ones(sim_len)      # Voltage input
time    = 1*np.arange(0, sim_len)   # Time-steps

# Counters
dtl = 0.0   # Time since last spike of left device
dtr = 0.0   # Time since last spike of right device

# Pre-factors
a   = (R0*(Rs*Cp + dt) + dt*Rs)/((Rs*Cp + dt)*R0)
b   = Rs*Cp/(Rs*Cp + dt)

# Lists
I       = []    # Total current
Il      = []    # Left branch current
Ir      = []    # Right branch current
Il_out  = []    # Left device output current
Ir_out  = []    # Right device output current
I0      = []    # Coupling capacitor current
Vcl     = []    # Left sample capacitor voltage
Vcr     = []    # Right sample capacitor voltage
V0      = []    # Coupling capacitor voltage
isil    = []    # Inter spike intervals of left device
isir    = []    # Inter spike intervals of right device

# Compute current of coupling capacitor   
def getI0():
    return 0. if C0 == 0. else (Vcr[-1] - Vcl[-1] - V0[-1])/(+ dt/C0 + 2*b*dt/(Cp*a))
    
# Compute voltage of coupling capacitor    
def getV0():
    return 0. if C0 == 0. else  V0[-1] + I0[-1]*dt/C0

# Compute voltage of left capacitor    
def getVcl(t):
    return b/a*(Vcl[-1] + (Vapp[t]/R0 + I0[-1])*dt/Cp)

# Compute voltage of right capacitor    
def getVcr(t):
    return b/a*(Vcr[-1] + (Vapp[t]/R0 - I0[-1])*dt/Cp)
    
# Compute current of left device (capacitor + sample)
def getIl(t):
    return (Vapp[t] - Vcl[-1])/R0 + I0[-1]  
 
# Compute current of right device (capacitor + sample) 
def getIr(t):
    return (Vapp[t] - Vcr[-1])/R0 - I0[-1]

# Compute probability of firing   
def getP(V, ts, V_th):
    if (V < V_th):
        return 1. - 1./(1. + ts/dt*np.exp((V - V_th)/dV))  
    else: # This helps with overflow errors
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
    Il_out.append(Il[-1] + I0[-1])
    Ir_out.append(Ir[-1] - I0[-1])
    
# Create 2D plot 
def makeFig(x, y, color, xLabel, yLabel, name, yLimBottom = None, yLimTop = None):    
    plt.figure()
    fig, ax = plt.subplots(figsize=(15, 4))
    ax.plot(x, y, color, marker='s')
    ax.set(xlabel=xLabel, ylabel=yLabel)
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)
    if yLimBottom is not None:
        ax.set_ylim(bottom=yLimBottom)
    if yLimTop is not None:
        ax.set_ylim(top=yLimTop)
    ax.tick_params(axis = 'both', which = 'both', labelsize = 16, length = 4)
    fig.tight_layout()
    plt.savefig('./figures/' + name + '_C0=' + str(C0) + '.pdf') 
    plt.close()

# Create histogram    
def makeHist(dist, width, color, xLabel, yLabel, name, 
            xLimLeft = None, xLimRight = None, yLimBottom = None, yLimTop = None):    
    plt.figure()
    fig, ax = plt.subplots(figsize=(15, 4)) 
    plt.hist(dist, bins=np.arange(min(dist), max(dist) + width, width), color=color, alpha=0.5)
    ax.set(xlabel=xLabel, ylabel=yLabel)
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)
    
    if yLimBottom is not None:
        ax.set_ylim(bottom=yLimBottom)
    if yLimTop is not None:
        ax.set_ylim(top=yLimTop)
    if xLimLeft is not None:
        ax.set_xlim(left=xLimLeft)
    if xLimRight is not None:
        ax.set_xlim(right=xLimRight)
        
    ax.tick_params(axis = 'both', which = 'both', labelsize = 16, length = 4)
    fig.tight_layout()
    plt.savefig('./figures/' + name + '_C0=' + str(C0) + '.pdf') 
    plt.close()

for C0 in C0Array:
    print('C0=' + str(C0))
    for t in time:   
        # Initial conditions
        if len(Vcl) == 0:
            Vcl.append(Vcl_0)
            Vcr.append(Vcr_0)
            V0.append(Vcl_0 - Vcr_0)
        
        solveCircuit(t)
        
        # Relaxation of the device to the insulating phase
        Vcl[-1] = Vcl[-1] if dtl > tRef else Vcl_0
        Vcr[-1] = Vcr[-1] if dtr > tRef else Vcr_0
                 
        # Spiking
        rand = r.random()    
        if rand < getP(Vcl[-1], dtl, Vl_th): 
            isil.append(dtl)
            dtl = 0.        
            Vcl[-1] = Vcl_0
        else:
            dtl = dtl + 1.
            
        if (rand < getP(Vcr[-1], dtr, Vr_th)):  
            isir.append(dtr)
            dtr = 0.        
            Vcr[-1] = Vcr_0
        else:
            dtr = dtr + 1.
            
        t = t + 1
    
    '''
    # Plot voltage
    makeFig(time, Vcl[:sim_len], 'black', 'Time (arb. units)', 'Voltage (arb. units)', 'Vcl', 0.0, 0.5)
    makeFig(time, Vcr[:sim_len], 'red', 'Time (arb. units)', 'Voltage (arb. units)', 'Vcr', 0.0, 0.5)
    
    # Plot current
    makeFig(time, Il_out, 'black', 'Time (arb. units)', 'Ouput current (arb. units)', 'Il_out', 0.0)
    makeFig(time, Ir_out, 'red', 'Time (arb. units)', 'Output current (arb. units)', 'Ir_out', 0.0)
    '''
    
    # Plot ISI histograms
    makeHist(isil, 100*dt, 'black', 'Inter spike intervals (arb. units)', 'Counts', 'ISI_l', 0.0, 150)
    makeHist(isir, 100*dt, 'red', 'Inter spike intervals (arb. units)', 'Counts', 'ISI_r', 0.0, 150)
    
    # Reset 
    dtl = 0.0   
    dtr = 0.0     
    I       = []    
    Il      = []   
    Ir      = []    
    Il_out  = []    
    Ir_out  = []    
    I0      = []    
    Vcl     = []    
    Vcr     = []    
    V0      = []    
    isil    = []    
    isir    = []  
    
    
    
    