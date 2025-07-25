import plot_tools as pt
import constants as c
import numpy as np
import matplotlib.pyplot as plt
import dynamics as dyn
from scipy import integrate as integ
np.set_printoptions(linewidth=175)

import os
sol_dir = os.path.join(os.getcwd(),"OrbitSolutions")


## Extract state information from loaded solution set
Sol = np.loadtxt('OrbitSolutions/Halo.csv', delimiter=',')
states = Sol[1:7,:]
tvec = Sol[0,:]
eps = 1*10**(-4)


Monodromy = Sol[7:,-1].reshape(6,6).transpose() # Monodromy matrix is STM after one full period of periodic orbit

eigvals, eigvecs = np.linalg.eig(Monodromy) # Calculate eigen-vectors and values of monodromy matrix
print(eigvals, '\n')
print(eigvecs, '\n')

###############
eig = 1 # Select the desired (correct) eigen value and vector to create perturbation vector
lambdaU = eigvals[eig] # Get the desired eigenvalue
vU = eigvecs[:,eig].reshape(6,1)              # Get the desired eigenvector for eigval of 1.026


## Create initial conditions of all the sampled states of the orbit
samples = 20 # Select number of trajectories to sample from orbit to create manifold    #round(np.size(states,axis=1)/200)
h = int(np.floor((len(tvec)-1)/samples)) # Iteration step value. Used to properly sample across periodic orbit
pertVecs = np.zeros((6, 2*samples)) # Pre-allocate matrix sized to hold each perturbed initial state vector. Each is 6x1

i=0
# Get delta vectors at a sample of points in the orbit
for point in np.arange(0,int(h*samples), h):   # at each time from original solved orbit, get the STM, multiply it by vU and scale, 
    STM = Sol[7:,point].reshape(6,6) # Extract STM at sample point and resize
    delta = np.matmul(STM,vU) # Delta vector is STM*vU at sampled point in orbit

    xUp = states[:,point].reshape(6,1) - eps* delta/np.linalg.norm(delta) # Create perturbed initial conditions at sampled point
    # xUn = states[:,point].reshape(6,1) - eps* delta/np.linalg.norm(delta)
    pertVecs[:,i] = xUp[:,0] # Save perturbed IC state vector to big matrix for convenience
    # pertVecs[:,i+1] = xUn[:,0]
    i+=1

# print(pertVecs)


## Integrate each perturbation vector in time. These with be the (un)stable trajectories from the periodic orbit.
# Setup the time bounds of each new trajectory
t0=0; T= 1.5*tvec[-1]
tsteps=1000
tspan=np.linspace(0, T,tsteps)
step = (T-t0)/tsteps


fig = plt.figure()
ax = plt.axes(projection='3d')
for point in range(samples): # Loop through each perturbed IC state vector
    state0syn = pertVecs[:,point]
    # print(state0syn)
    SynSol = integ.solve_ivp(fun=dyn.BackwardEOMs, t_span=[t0,T], y0=state0syn, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6) # Integrate the IC with the proper time direction EOMs.

    # Parse solution
    solvec = SynSol.y
    time = SynSol.t
    traj = pt.PlotManifold(solvec, time, c.mustar, ax, title='Stable Invariant Manifold', eigval=lambdaU, eigvec=vU) # Add propogated trajectory to plot
    # plt.show()

# ax.text2D(-.5, 0.5, (lambdaU, vU), transform=ax.transAxes)
# plt.show()

# pt.Orbit3D(states, tvec, c.mustar, args={'Frame':'Synodic'})                 # plot 3d orbit in synodic frame
# pt.PhasePortraits(states, tvec)
plt.show()


# sailAcc = np.zeros((3,len(SynSol.t)))
# for iter in range(len(SynSol.t)):
#     sailAcc[:,iter] = np.array([c.a_0*np.cos(c.ws*SynSol.t[iter]), -c.a_0*np.sin(c.ws*SynSol.t[iter]), 0])
# Jacobi = ot.JacobiConstant(SynSol.y, sailAcc)
# pt.JacobiPlot(Jacobi, SynSol.t, [-2.175,-2.05])
# print(ot.closeApproach(SynSol))
