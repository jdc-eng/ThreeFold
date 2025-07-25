"""CR3BP Dynamics File."""
import numpy as np
import constants as c

'''Sub-script of 1 refers to the larger primary, Earth. Sub-script 2 refers to the smaller primary, the Mooon.'''
def SynodicEOMs(t, state):
    x, y, z, vx, vy, vz = state[0:6]
    STM0 = state[6:].reshape(6,6)
    mu = c.mustar

    m1 = 1 - mu
    m2 = mu

    rb1 = x + mu
    rb2 = x + mu - 1

    r1 = np.sqrt( (rb1)**2 + y**2 + z**2 )          # Magnitude of r1-sat vector
    r2 = np.sqrt( (rb2)**2 + y**2 + z**2 )          # Magnitude of r2-sat vector

    statedot = np.zeros((42))

    statedot[0:3] = state[3:6]
    statedot[  3] = x + 2*vy - m1* (rb1)/ r1**3 - m2* (rb2)/ r2**3
    statedot[  4] = y - 2*vx - m1*     y/ r1**3 - m2*     y/ r2**3
    statedot[  5] =          - m1*     z/ r1**3 - m2*     z/ r2**3

    Uxx = 1 - (1-mu)*(r1**2-3*rb1**2)/r1**5 - mu*(r2**2-3*rb2**2)/r2**5
    Uxy = 3*y*(1-mu)*rb1/r1**5 + 3*mu*y*rb2/r2**5
    Uxz = 3*z*(1-mu)*rb1/r1**5 + 3*mu*z*rb2/r2**5
    Uyx = 3*y*(1-mu)*rb1/r1**5 + 3*mu*y*rb2/r2**5   # = Uxy
    Uyy = 1 - (1-mu)*(r1**2-3*y  **2)/r1**5 - mu*(r2**2-3*y  **2)/r2**5
    Uyz = 3*  (1-mu)*y*z/r1**5 + 3*mu*y*z/r2**5
    Uzx = 3*z*(1-mu)*rb1/r1**5 + 3*mu*z*rb2/r2**5   # = Uxz
    Uzy = 3*  (1-mu)*y*z/r1**5 + 3*mu*y*z/r2**5     # = Uyz
    Uzz =    -(1-mu)*(r1**2-3*z  **2)/r1**5 - mu*(r2**2-3*z  **2)/r2**5

    Z = np.zeros((3,3))
    I = np.eye(3,3)
    U = np.array([[Uxx, Uxy, Uxz],[Uyx, Uyy, Uyz],[Uzx, Uzy, Uzz]])
    Omega = np.array([[0, 2, 0],[-2, 0, 0],[0, 0, 0]])
    Amat = np.concatenate((np.concatenate((Z,I),axis=1),np.concatenate((U,Omega),axis=1)),axis=0)

    STM = np.matmul(Amat, STM0).reshape(1,36)
    statedot[6:] = STM

    return np.array(statedot)

def ForwardEOMs(t, state):
    x, y, z, vx, vy, vz = state
    mu = c.mustar

    m1 = 1 - mu
    m2 = mu

    rb1 = -mu
    rb2 = 1 - mu

    r1 = np.sqrt( (x - rb1)**2 + y**2 + z**2 )          # Magnitude of r1-sat vector
    r2 = np.sqrt( (x - rb2)**2 + y**2 + z**2 )          # Magnitude of r2-sat vector

    statedot = np.zeros(len(state))

    statedot[0:3] = state[3:6]
    statedot[3]   = x + 2*vy - m1* (x-rb1)/ r1**3 - m2* (x-rb2)/ r2**3
    statedot[4]   = y - 2*vx - m1*    y   / r1**3 - m2*    y   / r2**3
    statedot[5]   =           -m1*    z   / r1**3 - m2*    z   / r2**3


    return np.array(statedot)

def BackwardEOMs(t, state):
    x, y, z, vx, vy, vz = state
    mu = c.mustar

    m1 = 1 - mu
    m2 = mu

    rb1 = -mu
    rb2 = 1 - mu

    r1 = np.sqrt( (x - rb1)**2 + y**2 + z**2 )          # Magnitude of r1-sat vector
    r2 = np.sqrt( (x - rb2)**2 + y**2 + z**2 )          # Magnitude of r2-sat vector

    statedot = np.zeros(len(state))

    statedot[0:3] = -state[3:6]
    statedot[3]   = -(x + 2*vy - m1* (x-rb1)/ r1**3 - m2* (x-rb2)/ r2**3)
    statedot[4]   = -(y - 2*vx - m1*    y   / r1**3 - m2*    y   / r2**3)
    statedot[5]   = -(          -m1*    z   / r1**3 - m2*    z   / r2**3)


    return np.array(statedot)


def BarycentricEOMs(t, state):
    '''State vector is now 8x1; 7th element is theta=non-dim time'''
    x, y, z, vx, vy, vz = state
    mustar = c.mustar

    mu1 = 1 - mustar
    mu2 = mustar

    rb1 = - mustar
    rb2 = 1 - mustar

    r1 = np.sqrt( (x - rb1*np.cos(t))**2 + (y - rb1*np.sin(t))**2 + z**2 )      # Magnitude of r1-sat vector
    r2 = np.sqrt( (x - rb2*np.cos(t))**2 + (y - rb2*np.sin(t))**2 + z**2 )      # Magnitude of r2-sat vector

    statedot = np.zeros(len(state))
    statedot[0:3] = state[3:6]
    statedot[3] = -mu1*(x - rb1*np.cos(t))/r1**3 - mu2*(x + rb2*np.cos(t))/r2**3
    statedot[4] = -mu1*(y - rb1*np.sin(t))/r1**3 - mu2*(y - rb2*np.sin(t))/r2**3
    statedot[5] = -mu1*z/r1**3 - mu2*z/r2**3

    return np.array(statedot)


def SailSynodicEOMs(t, state):
    x, y, z, vx, vy, vz = state
    mustar = c.mustar
    mu4 = c.sunMass/(c.earthMass+c.moonMass)
    L = c.AU/c.lstar
    
    m1 = 1 - mustar
    m2 = mustar

    rb1 = -mustar
    rb2 = 1 - mustar

    r1 = np.sqrt( (x - rb1)**2 + y**2 + z**2 )          # Magnitude of r1-sat vector
    r2 = np.sqrt( (x - rb2)**2 + y**2 + z**2 )          # Magnitude of r2-sat vector


    # getJC = ot.JacobiConstant(state, np.zeros((3,1)))
    # if getJC < -2.96 : 
    acc_sail = np.array([c.a_0*np.cos(c.ws*t), -c.a_0*np.sin(c.ws*t), 0])
    # else:
    #     acc_sail = np.zeros((3,1))



    # a = L*np.cos(c.ws*t)
    # b = L*np.sin(c.ws*t)
    # acc_sun = np.array([-mu4*a/ (a**2 + b**2 )**(3/2) - mu4*(x-a)/( (x-a)**2 + (y+b)**2 + z**2)**(3/2),
    #                    mu4*b/ (a**2 + b**2 )**(3/2) - mu4*(y+b)/( (x-a)**2 + (y+b)**2 + z**2)**(3/2),
    #                   -mu4*z/( (x-a)**2 + (y+b)**2 + z**2)**(3/2)])
    

    statedot = np.zeros(len(state))
    statedot[0:3] = state[3:]
    statedot[3] = x + 2*vy - m1* (x-rb1)/ r1**3 - m2* (x-rb2)/ r2**3 + acc_sail[0] # + acc_sun[0]
    statedot[4] = y - 2*vx - m1*    y   / r1**3 - m2*    y   / r2**3 + acc_sail[1] # + acc_sun[1]
    statedot[5] =           -m1*    z   / r1**3 - m2*    z   / r2**3 + acc_sail[2] # + acc_sun[2]

    return np.array(statedot)


def StateTransMat(t, state, STM0):
    '''state is '''
    mu = c.mustar
    x, y, z = state[0:3]
    
    Uxx = -((-2*(x + mu)**2 + y**2 +   z**2)*(1 - mu)) / ((x + mu)**2 + y**2 + z**2)**(5/2) - (mu*(y**2 + z**2 - 2*(x + mu - 1)**2)) / ((x + mu - 1)**2 + y**2 + z**2)**(5/2) + 1
    Uxy = (3*y*(-mu + 1)*(mu + x)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) + (3*mu*y*(mu + x - 1)) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))
    Uxz = (3*z*(-mu + 1)*(mu + x)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) + (3*mu*z*(mu + x - 1)) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))
    Uyx = (3*y*(-mu + 1)*(mu + x)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) + (3*mu*y*(mu + x - 1)) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))
    Uyy = -((   (x + mu)**2 + z**2 - 2*y**2)*(1 - mu)) / ((mu + x)**2 + y**2 + z**2)**(5/2) - (mu*(-2*y**2 + mu**2 - 2*mu + x**2 + 2*mu*x - 2*x + z**2 + 1)) / ((mu + x - 1)**2 + y**2 + z**2)**(5/2) + 1
    Uyz = (3*y*z*(-mu + 1)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) + (3*mu*y*z) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))
    Uzx = (3*z*(-mu + 1)*(x + mu)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) + (3*mu*z*(x + mu - 1)) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))
    Uzy = (3*y*z*(-mu + 1)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) + (3*mu*y*z) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))
    Uzz = -((   (x + mu)**2 + y**2 - 2*z**2)*(1 - mu)) / ((y**2 + z**2 + (mu + x)**2)**(5/2)) -(mu*(-2*z**2 + mu**2 - 2*mu + y**2 + x**2 + 2*mu*x - 2*x + 1)) / ((y**2 + z**2 + (mu + x - 1)**2)**(5/2))

    Z = np.zeros((3,3))
    I = np.eye(3,3)
    U = np.array([[Uxx, Uxy, Uxz],[Uyx, Uyy, Uyz],[Uzx, Uzy, Uzz]])
    Omega = np.array([[0, 2, 0],[-2, 0, 0],[0, 0, 0]])

    Amat = np.concatenate((np.concatenate((Z,I),axis=1),np.concatenate((U,Omega),axis=1)),axis=0)

    STM = np.matmul(Amat, STM0)


    return STM



def two_body_ode( t, state, mu = c.EMmu ):
	# state = [ rx, ry, rz, vx, vy, vz ]

	r = state[ :3 ]
	a = -mu * r / np.linalg.norm( r ) ** 3

	return np.array( [
		state[ 3 ], state[ 4 ], state[ 5 ],
		    a[ 0 ],     a[ 1 ],     a[ 2 ] ] )