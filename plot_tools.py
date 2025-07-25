"""Plotting tools file."""
import numpy as np
import matplotlib.pyplot as plt
import constants as c

## honestly should replace most of this with a plotting class that lets me do shit with a lot more abstraction

def genPlot(xvals, yvals, xbound, ybound, labels={'xlabel':'Time (s)','ylabel':'Position (km)'}, legend={0:'Label1'}):
    # numplot = np.shape(yvals)[1]
    # print(numplot)
    # axs = plt.plot() 
    fig1, axs = plt.subplots()
    # i = 0
    # for signal in yvals[:,0]:
    plt.plot(xvals, yvals)
        # i += 1

    # plt.plot(xvals, yvals[1,:])
    # axs[0].axis('equal')
    axs.set_title('X-Pos v. Time (s)', fontsize=10)
    axs.set_xlim(xbound[0], xbound[1])
    axs.set_ylim(ybound[0], ybound[1])
    plt.xlabel(labels['xlabel'])
    plt.ylabel(labels['ylabel'])


def compPlot(tvec, solvec):
    ''' tvec is odesol.t
        solvec is the odesol.y'''

    x = np.array(solvec[0,:])
    y = np.array(solvec[1,:])
    z = np.array(solvec[2,:])

    fig1, axs = plt.subplots(3,1)

    axs[0].plot(tvec, x)
    # axs[0].axis('equal')
    axs[0].set_title('X-Pos v. Time (s)', fontsize=10)
    axs[0].set_xlim(tvec[0], tvec[-1])
    # axs[0].set_ylim(np.min(x), np.max(x))

    axs[1].plot(tvec, y)
    # axs[1].axis('equal')
    axs[1].set_title('Y-Pos v. Time (s)', fontsize=10)
    axs[1].set_xlim(tvec[0], tvec[-1])
    axs[1].set_ylim(np.min(y), np.max(y))

    axs[2].plot(tvec, z)
    # axs[2].axis('equal')
    axs[2].set_title('Z-Pos v. Time (s)', fontsize=10)
    axs[2].set_xlim(tvec[0], tvec[-1])
    # axs[2].set_ylim(np.min(z), np.max(z))

    # return fig1

def PhasePortraits(solvec, tvec):

    x    = np.array(solvec[0,:])
    y    = np.array(solvec[1,:])
    z    = np.array(solvec[2,:])
    xdot = np.array(solvec[3,:])
    ydot = np.array(solvec[4,:])
    zdot = np.array(solvec[5,:])

    fig = plt.figure(2)
    ax1 = plt.axes()
    traj = ax1.scatter(x,xdot, c=tvec, cmap = 'plasma', s=0.5)
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Xdot')
    plt.colorbar(traj)

    plt.figure(3)
    ax2 = plt.axes(projection='3d')
    traj2 = ax2.scatter(y,ydot,x, c=tvec, cmap = 'plasma', s=0.5)
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Xdot')
    plt.colorbar(traj2)
    


def Orbit3D(solvec, time, mu, args={}):
    _args = {'Frame': 'Synodic'}
    x = np.array(solvec[0,:])
    y = np.array(solvec[1,:])
    z = np.array(solvec[2,:])

    
    ax = plt.axes(projection='3d')
    traj = ax.scatter(x,y,z, c=time, cmap = 'plasma', s=.5)
    ax.scatter(0,0,0, c='m', marker='*')

    if args['Frame'] == 'ECRF':
        n = np.linspace(0,2*np.pi,len(time))
        v = np.linspace(0, np.pi, len(time))

        re = c.earthD / c.lstar
        rm = c.moonD / c.lstar

        EarthX = -mu*np.cos(n)
        EarthY = -mu*np.sin(n)
        MoonX = (1)*np.cos(n)
        MoonY = (1)*np.sin(n)

        xe = re * np.outer(np.cos(n), np.sin(v)) 
        ye = re * np.outer(np.sin(n), np.sin(v)) + 0
        ze = re * np.outer(np.ones(np.size(n)), np.cos(v)) + 0

        xm = rm * np.outer(np.cos(n), np.sin(v)) + (1)
        ym = rm * np.outer(np.sin(n), np.sin(v))
        zm = rm * np.outer(np.ones(np.size(n)), np.cos(v))
        # c=time, cmap = 'plasma',
        ax.scatter(EarthX, EarthY, 0, c='g', s=1)
        ax.scatter(MoonX, MoonY, 0, c='m', s=1)
        ax.plot_surface(xe,ye,ze)
        ax.plot_surface(xm,ym,zm)
        plt.title('Orbit in the Earth Centered Inertial Frame (ECRF).')
    
    else: 
        n = np.linspace(0,2*np.pi,100)
        v = np.linspace(0, np.pi, 100)

        re = c.earthD / c.lstar;
        rm = c.moonD / c.lstar;

        xe = re * np.outer(np.cos(n), np.sin(v)) - mu
        ye = re * np.outer(np.sin(n), np.sin(v)) + 0
        ze = re * np.outer(np.ones(np.size(n)), np.cos(v)) + 0

        xm = rm * np.outer(np.cos(n), np.sin(v)) + (1-mu)
        ym = rm * np.outer(np.sin(n), np.sin(v))
        zm = rm * np.outer(np.ones(np.size(n)), np.cos(v))

        ax.plot_surface(xe,ye,ze)
        ax.plot_surface(xm,ym,zm)
        plt.title('Orbit in the Earth-Moon Rotating Frame')

    plt.axis('equal')
    ax.legend()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar(traj)


def PlotManifold(solvec, time, mu, ax, title,eigval, eigvec):
    # _args = {'Frame': 'Synodic'}
    x_vals = np.array(solvec[0,:])
    y_vals = np.array(solvec[1,:])
    z_vals = np.array(solvec[2,:])

    
    # ax = plt.axes(projection='3d')
    traj = ax.scatter(x_vals,y_vals,z_vals, c=time, cmap = 'plasma',s=.5)
    ax.scatter(0,0,0, c='m', marker='*')

    n = np.linspace(0,2*np.pi,100)
    v = np.linspace(0, np.pi, 100)

    re = c.earthD / c.lstar
    rm = c.moonD / c.lstar

    xe = re * np.outer(np.cos(n), np.sin(v)) - mu
    ye = re * np.outer(np.sin(n), np.sin(v)) + 0
    ze = re * np.outer(np.ones(np.size(n)), np.cos(v)) + 0

    xm = rm * np.outer(np.cos(n), np.sin(v)) + (1-mu)
    ym = rm * np.outer(np.sin(n), np.sin(v))
    zm = rm * np.outer(np.ones(np.size(n)), np.cos(v))

    ax.plot_surface(xe,ye,ze)
    ax.plot_surface(xm,ym,zm)
    plt.suptitle(title)

    ax.set_title(("Eigenvalue: ", eigval ))
    plt.axis('equal')
    # ax.text2D(0.05, 0.95, (r'Eigenvalue: ', eigval, r'\nEigvec: ', eigvec), transform=ax.transAxes)
    ax.legend()
    plt.xlabel('X\n')
    plt.ylabel('Y\n')
    # plt.colorbar(traj)


def Orbit2D(solvec, time, mu, args={}):
    _args = {'Frame':'Synodic'}
    x = np.array(solvec[0,:])
    y = np.array(solvec[1,:])
    z = np.array(solvec[2,:])

    fig = plt.figure()
    ax = plt.axes()
    traj = ax.scatter(x,y, c=time, cmap = 'plasma', s=.5, label='Spacecraft')
    ax.plot(-mu, 0, c='g', marker='x', label='Earth')
    ax.plot(1-mu, 0, c='b', marker='^', label='Moon')
    ax.plot(0,0, c='m', marker='*')


    # if args['Frame'] == 'Barycentric':
    #     n = np.linspace(0,2*np.pi,100)
    #     v = np.linspace(0, np.pi, 100)

    #     re = c.earthD / c.lstar
    #     rm = c.moonD / c.lstar

    #     EarthX = -mu*np.cos(n)
    #     EarthY = -mu*np.sin(n)
    #     MoonX = (1-mu)*np.cos(n)
    #     MoonY = (1-mu)*np.sin(n)

    #     xe = re * np.outer(np.cos(n), np.sin(v)) + mu
    #     ye = re * np.outer(np.sin(n), np.sin(v)) + 0

    #     xm = rm * np.outer(np.cos(n), np.sin(v)) - (1-mu)
    #     ym = rm * np.outer(np.sin(n), np.sin(v))

    #     ax.scatter(EarthX, EarthY, c='g')
    #     ax.scatter(MoonX, MoonY, c='b')
    #     ax.scatter(xe,ye)
    #     ax.scatter(xm,ym, c=time, cmap = 'plasma', label='Moon Orbit')

    plt.axis('equal')
    ax.legend()
    ax.set_xlim(-1,1)
    plt.xlabel('X')
    plt.ylabel('Y')
    # plt.colorbar(traj)


def JacobiPlot(C, time, ybound ):
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(time, C)

    # ax.set_xlim(xbound[0], xbound[1])
    ax.set_ylim(ybound[0], ybound[1])
    plt.title('Jacobi Constant vs. Time')
    plt.xlabel('Time')
    plt.ylabel('Jacobi Constant C')


def PhasePortrait(solvec, tvec):
    ''' tvec is odesol.t
        solvec is the odesol.y'''

    x,y,z,vx,vy,vz = solvec

    fig1, axs = plt.subplots(3,1)

    axs[0].plot(x, vx)
    axs[0].plot(x[0], vx[0], '*','tab:orange')
    axs[0].set_title('X-Vel v. X-Pos', fontsize=10)
    # axs[0].set_xlim(x[0], x[-1])
    # axs[0].set_ylim(np.min(x), np.max(x))
    # axs[0].imshow(c=tvec, cmap='plasma')

    axs[1].plot(y, vy)
    axs[1].plot(y[0], vy[0], '*','tab:orange')
    # axs[1].axis('equal')
    axs[1].set_title('Y-Vel v. Y-Pos', fontsize=10)
    # axs[1].set_xlim(y[0], y[-1])
    # axs[1].set_ylim(np.min(y), np.max(y))

    axs[2].plot(z, vz)
    axs[2].plot(z[0], vz[0], '*','tab:orange')
    # axs[2].axis('equal')
    axs[2].set_title('Z-Vel v. Z-Pos', fontsize=10)
    # axs[2].set_xlim(z[0], z[-1])
    # axs[2].set_ylim(np.min(z), np.max(z))

    # return fig1

