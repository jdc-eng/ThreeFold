"""Orbit calculations and reference frames tools. Mostly functions for frame transformations."""
'''Thanks Niko, I bow to you.'''
import constants as c
import numpy as np
import dynamics
from scipy import integrate as integ



def invertKeplerTimeEquation(M,E0,e):
    """Invert the Kepler Time Equation
    INPUTS:
        M - mean anomaly
        E0 - initial guess for Eccentric Anomaly
        e - eccentricity
    OUTPUTS:
        E - Eccentric Anomaly
    """
    E=E0
    while (np.abs(M-E+e*np.sin(E))>1e-12):
        E = E - (M-E+e*np.sin(E))/(e*np.cos(E)-1)
    return E

def kepler2rv(a,e,Omega,I,omega,tmtp,mu):
    """ Convert Keplerian Elements to a position and velocity.
    INPUTS:
        a - semjax 
        e - eccentricity
        Omega - longitude of ascending node
        I - inclination
        omega - argument of periapsis
        tmtp - time since periapsis passage (t-tp)
        mu - gravitational parameter (units of this dictate units of output)
    OUTPUTS:
        r - inertial position vector
        v - inertial velocity vector
    """
    b = a*np.sqrt(1-e**2)
    Omega_rotation = np.array([[np.cos(Omega), np.sin(Omega), 0],[-np.sin(Omega), np.cos(Omega), 0],[0, 0, 1]])
    I_rotation = np.array([[1, 0, 0],[0, np.cos(I), np.sin(I)],[0, -np.sin(I), np.cos(I)]])
    omega_rotation = np.array([[np.cos(omega), np.sin(omega), 0],[-np.sin(omega), np.cos(omega), 0],[0, 0, 1]])
    pDCMi = np.matmul(omega_rotation,np.matmul(I_rotation,Omega_rotation))
    n = np.sqrt(mu/a**3)
    M = np.mod(n*(tmtp),2*np.pi)
    E = invertKeplerTimeEquation(M,M,e)
    r_p = np.transpose(np.array([a*(np.cos(E)-e),b*np.sin(E),0]))
    v_p = (a*n/np.linalg.norm(r_p))*np.array([-a*np.sin(E),b*np.cos(E),0])
    r = np.matmul(pDCMi.transpose(),r_p)
    v = np.matmul(pDCMi.transpose(),v_p)
    return r,v

def rv2kepler(r,v,mu):
    """ Convert a position and velocity to Keplerian Elements.
    INPUTS:
        r - inertial position vector
        v - inertial velocity vector
        mu - gravitational parameter (units of this dictate units of output)
    OUTPUTS:
        a - semjax 
        e - eccentricity
        I - Inclination
    """
    rhat = r/np.linalg.norm(r)
    vhat = v/np.linalg.norm(v)
    h = np.cross(r,v)
    hhat = h/np.linalg.norm(h)

    ev = np.cross(v,np.cross(r,v))/mu-r/np.linalg.norm(r)
    e = np.linalg.norm(ev)
    ehat = ev/e

    nu = np.arctan2(np.dot(np.cross(ehat, rhat),hhat),np.dot(rhat, ehat))

    I = np.arccos(np.dot(hhat, np.array([0,0,1])))
    if e<1:
        a = np.linalg.norm(h)**2/(mu*(1-e**2))
    elif np.abs(e-1)< 1e-8:
        a = np.Infinity
    else:
        a = -np.linalg.norm(r)*(1+e*np.cos(nu))/(e**2-1)
    return a,e,I

def SolTransform(sol, rotation):
    _args = {'Bar2ECI': Bar2ECI,
             'Syn2Bar': Syn2Bar,
             'Bar2Syn': Bar2Syn,
             'ECI2Bar': ECI2Bar
            }
    rotState = np.zeros( (6, len(sol['t'])) )

    fun = _args[rotation]
    i = 0
    for t in sol['t']:
        rotState[0,i], rotState[1,i], rotState[2,i], rotState[3,i], rotState[4,i], rotState[5,i]= fun(sol['y'][:,i], t)
        i += 1
    RotSol = {'y':rotState, 't':sol['t']}
    return RotSol

def Syn2Bar(synstate, t, args={}):
    '''Assumes both states coincide at t=0'''
    x,y,z,vx,vy,vz = synstate
    rSYN = np.matrix([[x],[y],[z]])
    vSYN = np.matrix([[vx],[vy],[vz]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    J = np.matrix([[0,  1, 0],
                   [-1, 0, 0],
                   [0,  0, 0]])

    rBAR = np.matmul( At, rSYN)
    vBAR = -np.matmul( At, np.matmul( J, rSYN ) ) + np.matmul( At, vSYN)

    barstate = np.zeros(len(synstate))
    barstate[0], barstate[1], barstate[2] = rBAR
    barstate[3], barstate[4], barstate[5] = vBAR
    return barstate

def Bar2Syn(barstate, t):
    x,y,z,vx,vy,vz = barstate
    rBAR = np.matrix([[x], [y], [z]])
    vBAR = np.matrix([[vx],[vy],[vz]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    
    J = np.matrix([[0,  1, 0],
                   [-1, 0, 0],
                   [0,  0, 0]])
    
    rEMR = np.matmul( np.transpose(At), rBAR )
    vEMR = np.matmul( np.transpose(At), ( vBAR + np.matmul( At, np.matmul( J, rEMR )) ) )

    synstate = np.zeros(len(barstate))
    synstate[0], synstate[1], synstate[2] = rEMR
    synstate[3], synstate[4], synstate[5] = vEMR
    return synstate

def Bar2ECI(state, t):
    mu = c.mustar
    rBAR = np.matrix([[state[0]], [state[1]], [state[2]]])
    vBAR = np.matrix([[state[3]], [state[4]], [state[5]]])

    rE   = np.matrix([[-c.mustar],       [0],        [0]])
    vE   = np.matrix([       [0],    [-mu*t],        [0]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    
    rECI = rBAR + np.matmul( At, rE )
    vECI = vBAR + np.matmul( At, vE )

    ECIstate = np.zeros(len(state))
    ECIstate[0], ECIstate[1], ECIstate[2] = rECI
    ECIstate[3], ECIstate[4], ECIstate[5] = vECI
    return ECIstate 

def ECI2Bar(state, t):
    mu = c.mustar
    rECI = np.matrix([[state[0]], [state[1]], [state[2]]])
    vECI = np.matrix([[state[3]], [state[4]], [state[5]]])

    rE   = np.matrix([[-mu],       [0],        [0]])
    vE   = np.matrix([       [0],    [-mu*t],        [0]])

    At = np.matrix([[np.cos(t), -np.sin(t), 0],
                    [np.sin(t),  np.cos(t), 0],
                    [        0,          0, 1]])
    
    rBAR = rECI - np.matmul( At, rE )
    vBAR = vECI - np.matmul( At, vE )

    BARstate = np.zeros(len(state))
    BARstate[0], BARstate[1], BARstate[2] = rBAR
    BARstate[3], BARstate[4], BARstate[5] = vBAR
    return BARstate



def computeManifold(Sol, samples, epsilon, T_factor, tsteps, args={}):
    '''Compute the appropriate manifold and return the 2d array of trajectory solutions. 
    Inputs: solution array, monodromy, stability, num_samples, epsilon, t0, t1, tsteps.  
    Sol: solution array of tvec, position, velocity, and STMs
    T_factor: scaling factor applied to final value in tvec. Used for time propogation extent of sampled ICs. Effectively ks number of periodic orbital periods to propogate each sampled state.
    tsteps: number of integration steps to take. Controls resolution of propogations.
    '''
    _args = {'Stability':'Stable'} # Default arg setting
    settings = {'Stable':{ 'index':1, 'EOMs':dynamics.BackwardEOMs, 'eps_switch':-1}, 
                'Unstable':{ 'index':0, 'EOMs':dynamics.ForwardEOMs, 'eps_switch':1} }
    
    index = settings[args['Stability']]['index']
    EOMs = settings[args['Stability']]['EOMs']
    eps_switch = settings[args['Stability']]['eps_switch']
    
    # Extract from Sol
    states = Sol[1:7,:]
    tvec = Sol[0,:]
    Monodromy = Sol[7:,-1].reshape(6,6).transpose()  # Monodromy matrix is STM after one orbital period


    eigvals, eigvecs = np.linalg.eig(Monodromy)      # Calculate eigen-vectors and values of monodromy matrix
    
    ## Collect proper eigens
    lambdaU = eigvals[index]                           # Get the desired eigenvalue
    vU = eigvecs[:,index].reshape(6,1)                 # Get the desired eigenvector for eigval of 1.026


    ## Get delta vectors for initial conditions of sampled states
    h = int(np.floor((len(tvec)-1)/samples))         # Orbital state step value. Used to properly sample across periodic orbit
    pertVecs = np.zeros((6, samples))              # Pre-allocate matrix sized to hold each perturbed initial state vector. Each is 6x1

    i=0
    for point in np.arange(0,int(h*samples), h):     # at each time from original solved orbit, get the STM, multiply it by vU and scale,
        STM = Sol[7:,point].reshape(6,6)             # Extract STM at sample point and resize
        delta = np.matmul(STM,vU)                    # Delta vector is STM*vU at sampled point in orbit

        xU = states[:,point].reshape(6,1) +  epsilon* delta/np.linalg.norm(delta) # Create perturbed initial conditions at sampled point
        pertVecs[:,i] = xU[:,0]                     # Save perturbed IC state vector to big matrix for convenience
        # pertVecs[:,i+1] = xUn[:,0]
        i+=1
    
    
    ## Integrate pert vectors
    t0 = 0; tf = T_factor*tvec[-1]
    tspan=np.linspace(t0, tf, tsteps)
    step = (tf-t0)/tsteps
    
    # Pre-create 3D array
    # states, steps = np.shape(Sol)
    manifold = np.zeros((samples, 7, tsteps))
    for point in range(samples): # Loop through each perturbed IC state vector
        state0 = pertVecs[:,point]
        # print(state0syn)
        SynSol = integ.solve_ivp(fun=EOMs, t_span=[t0,tf], t_eval=tspan, y0=state0, method='DOP853', max_step = step, atol=1e-9, rtol=1e-6) # Integrate the IC with the proper time direction EOMs.
        manifold[point, 0,:] = SynSol.t
        manifold[point,1:,:] = SynSol.y[0:6,:]
            
    

    return manifold




def Tisserand(state,mu):
    ## Convert to regular shit first
    r = state.y[0:3]* c.lstar
    v = state.y[3:]* (c.lstar/c.tstar)
    # print(r)
    TissCrit = np.zeros(len(state.t))
    
    for i in range(len(state.t)) :
        # print(r[:,i]) 
        a, e, I = rv2kepler(r[:,i],v[:,i],mu)
        
        TissCrit[i] = 1/a + a*np.sqrt(a*(1-e**2))*np.cos(I)
    return TissCrit

def JacobiConstant(state, acc):
    '''Finds the Jacobi Constant for either one set of initial conditions or at each point in an orbit.'''
    
    mu = c.mustar
    rb1 = -mu
    rb2 = 1 - mu

    C = np.zeros(state.shape[1])
    for iter in range(np.shape(state)[1]):
        x,y,z,vx,vy,vz = state[:,iter]
        ax, ay, az = acc[:,iter]
        r1 = np.sqrt( (x - rb1)**2 + y**2 + z**2 )          # Magnitude of r1-sat vector
        r2 = np.sqrt( (x - rb2)**2 + y**2 + z**2 )          # Magnitude of r2-sat vector

        C[iter] = 2*(-.5*( x**2 + y**2 ) - ( (1-mu)/r1 + mu/r2 ) + .5*( vx**2 + vy**2 + vz**2 ) + ( ax*x + ay*y + az*z ))

    return C


def Dimensioning(pos, vel, acc, time):
    '''Returns the dimensioned or nondimensionalized values for the 3-body problem. '''


    return


def SailCharacteristics(a0):
    '''Returns the non-dimensionalized characteristic acceleration of the solar sail from mm/s^2'''

    a_0 = a0*10**(-6)/ ( c.lstar/c.tstar**2 )
    return a_0

def closeApproach(sol):
    """Returns the magnitude of the closest approach throughout the solution."""
    pos = sol.y[0:3]
    pmags = np.linalg.norm(pos, axis=0)
    posclose = np.min(pmags)

    return posclose, posclose*c.lstar

# def Syn2ECI(synstate, t):
#     '''Assumes both states coincide at t=0'''
#     x,y,z,vx,vy,vz = synstate
#     At = np.matrix([[np.cos(t), -np.sin(t), 0],
#                     [np.sin(t),  np.cos(t), 0],
#                     [        0,          0, 1]])

#     Vt = np.matrix([[vx - y],
#                     [vy + x + c.mustar],
#                     [  vz  ]])

#     r = np.matrix([[x+c.mustar],
#                    [y],
#                    [z]])

#     barstate = np.zeros(len(synstate))
#     barstate[0:3] = np.matmul( At, synstate[0:3] )
#     barstate[3], barstate[4], barstate[5]  = np.matmul( At, Vt)


# def ECI2Syn(state, t):
#     rECI = np.matrix([[state[0]], [state[1]], [state[2]]])
#     vECI = np.matrix([[state[3]], [state[4]], [state[5]]])
#     rE   = np.matrix([[c.mustar],[0],[0]])

#     At = np.matrix([[np.cos(t), -np.sin(t), 0],
#                     [np.sin(t),  np.cos(t), 0],
#                     [        0,          0, 1]])
    
#     J = np.matrix([[0,  1, 0],
#                    [-1, 0, 0],
#                    [0,  0, 0]])
    
#     rEMR =  np.matmul( np.transpose(At), rECI ) + rE
#     vEMR = -np.matmul( np.transpose(At), vECI ) - np.matmul(J, ( rEMR + rE ))
#     synstate = np.zeros(len(state))
#     synstate[0], synstate[1], synstate[2] = rEMR
#     synstate[3], synstate[4], synstate[5] = vEMR
#     return synstate


