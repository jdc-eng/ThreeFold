'''This file houses the simulation.'''
import constants as c
import orbit_tools as ot
import plot_tools as pt
import dynamics as cr3bp
from scipy import integrate as int
import numpy as np
import matplotlib.pyplot as plt
import state_library as slib
np.set_printoptions(linewidth=175)
## Setup timestep of propagation
t0 = 0                              # Initial tTime
tbound =  2.7764143751594395        # Final time
tsteps = 6000
step = (tbound-t0)/tsteps

## Create initial state vector
STM0 = np.eye(6,6).reshape(1,36)
state0syn = np.concatenate(([np.array(slib.StateDict('Halo'))], STM0), axis=1)[0]


## Integrate initial state with a certian integrator and dynamics model/reference frame
SynSol = int.solve_ivp(fun=cr3bp.SynodicEOMs, t_span=[t0,tbound], y0=state0syn, method='DOP853', max_step = step, atol=1e-12, rtol=1e-9)
# sailSol = int.solve_ivp(fun=cr3bp.SailSynodicEOMs, t_span=[t0,tbound], y0=state0syn, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6) 


Mono = SynSol.y[6:,-1].reshape(6,6)
print(Mono)
## Save to CSV file
output = np.concatenate((np.asmatrix(SynSol.t),SynSol.y), axis=0)
np.savetxt("OrbitSolutions/Halo.csv", output, delimiter=",")


## Plotting
pt.Orbit3D(SynSol.y, SynSol.t, c.mustar, args={'Frame':'Synodic'})                 # plot 3d orbit in synodic frame
plt.show()
