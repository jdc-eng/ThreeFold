import plot_tools
import orbit_tools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.widgets import Button, TextBox
import warnings
from matplotlib import cm


import os
sol_dir = os.path.join(os.getcwd(),"OrbitSolutions")

## Extract state information from loaded solution set
orbit = 'L1-Lyapunov.csv' # Select orbit solution file
global orbit_sol
orbit_sol = np.loadtxt(os.path.join(sol_dir, orbit), delimiter=',')
Monodromy = orbit_sol[7:,-1].reshape(6,6).transpose()  # Monodromy matrix is STM after one orbital period
eigvals, eigvecs = np.linalg.eig(Monodromy)      # Calculate eigen-vectors and values of monodromy matrix

lam1 = eigvals[0]; lam2 = eigvals[1]

#states = Sol[1:7,:]
#tvec = Sol[0,:]



## Setup manifold computation settings. Everything that has  text box needs to be set
# copy from SimPlots.y
global epsilon, samples, T_factor, tsteps
samples = 40
epsilon = 1*10**(-6)
T_factor = 1.5
tsteps = 2000



## Compute Monodromy atrix and eigen vals and vectors

# Add this near the top of the file, after imports but before any other functions
def add_gradient_background(ax):
    """Add a subtle gradient background to plots"""
    nx = 100
    ny = 100
    gradient = np.zeros((ny, nx, 4))
    gradient[:, :, 3] = np.linspace(0.2, 0, nx)
    ax.imshow(
        gradient,
        extent=[ax.get_xlim()[0], ax.get_xlim()[1], ax.get_ylim()[0], ax.get_ylim()[1]],
        aspect="auto",
        zorder=0,
    )


# Use a default matplotlib style that's always available
plt.style.use("default")

# Set some default plot parameters for better visualization
plt.rcParams["figure.figsize"] = [12, 8]
plt.rcParams["axes.grid"] = True
plt.rcParams["font.size"] = 12
plt.rcParams["lines.linewidth"] = 2

# Add these settings after the existing rcParams
plt.style.use("dark_background")
plt.rcParams.update(
    {
        "figure.facecolor": "#1e1e1e",
        "axes.facecolor": "#2d2d2d",
        "axes.edgecolor": "#666666",
        "axes.labelcolor": "white",
        "text.color": "white",
        "xtick.color": "white",
        "ytick.color": "white",
        "grid.color": "#444444",
        "legend.facecolor": "#2d2d2d",
        "legend.edgecolor": "#666666",
    }
)

# Define a modern color palette with more distinct colors
colors = ["#4cc9f0", "#f72585", "#ffd60a", "#7209b7", "#00f5d4", "#ff9e00"]
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=colors)


def create_dashboard(fig):
    fig.clear()
    global buttons 
    buttons = []
    
    global info_buttons
    info_buttons = []
    def create_orbit_tab(event=None):
        fig.clear()
        
        global buttons
        global info_buttons
        
        info_buttons.clear()

        def save_epsilon(epsilon_value):
            # set man_recompute to False initially?
            global epsilon
            print(epsilon)
            if float(epsilon_value) == epsilon:
                pass
            else:
                epsilon = float(epsilon_value)
            print(epsilon)
            return

        def save_samples(samples_value):
            # set man_recompute to False initially?
            global samples
            print(samples)
            if int(samples_value) == samples:
                pass
            else:
                samples = int(samples_value)
            print(samples)
            return
        
        def save_tfactor(tfactor_value):
            # set man_recompute to False initially?
            global T_factor
            print(T_factor)
            if float(tfactor_value) == T_factor:
                pass
            else:
                T_factor = float(tfactor_value)
            print(T_factor)
            return
        
        def save_tsteps(tsteps_value):
            # set man_recompute to False initially?
            global tsteps
            if int(tsteps_value) == tsteps:
                pass
            else:
                tsteps = int(tsteps_value)
            print(f'Changes number of timesteps to {tsteps}')
            return


        def plot_stable_manifold(label):
            '''Inputs: 3d trajectory solutions vector, plot title, eig value'''
            
            # clear subplot
            orbit_ax.clear()
            
            # plot stable manifold
            plot_tools.PlotManifold(orbit_ax, stable_manifold, line_color="green" )
            
            # re-add orbit
            plot_tools.Orbit3D(orbit_sol, orbit_ax, args={'Frame':'Synodic'})
            
            return
        
        def plot_unstable_manifold(label):
            '''Inputs: 3d trajectory solutions vector, plot title, eig value'''
            
            # clear subplot
            orbit_ax.clear()
            
            # plot stable manifold
            plot_tools.PlotManifold(orbit_ax, unstable_manifold, line_color="red" )

            # re-add orbit
            plot_tools.Orbit3D(orbit_sol, orbit_ax, args={'Frame':'Synodic'})
            
            return

        def refresh_manifolds(label):
            '''Recompute the stable and unstable manifold'''
            print("Computing manifolds.")
            global stable_manifold, unstable_manifold

            stable_manifold = orbit_tools.computeManifold(orbit_sol, samples, epsilon, T_factor, tsteps, args={'Stability':'Stable'})
            unstable_manifold = orbit_tools.computeManifold(orbit_sol, samples, epsilon, T_factor, tsteps, args={'Stability':'Unstable'})

            print("Completed manifold computation.")
            
            twod_plot.clear()
            plot_tools.PlotManifold2D(twod_plot, stable_manifold, line_color="green")
            plot_tools.PlotManifold2D(twod_plot, unstable_manifold, line_color="red")
            
            plot_tools.Orbit2D(twod_plot, orbit_sol)
            return
        

        gs = fig.add_gridspec(4, 2, top=0.9)
        
        ## Normal 3D orbit plot
        
        orbit_ax = fig.add_subplot(gs[0:3, 0], projection='3d')
        plot_tools.Orbit3D(orbit_sol, orbit_ax, args={'Frame':'Synodic'})

        twod_plot = fig.add_subplot(gs[0:3, 1])
        plot_tools.Orbit2D(twod_plot, orbit_sol)
        # twod_plot.set(xlim=(0.75,1.25))

        # subgrid for the buttons
        gs_info = gs[6].subgridspec(5, 2)
        
        # first eigen value
        eig_ax1= fig.add_subplot(gs_info[0,0])
        textstr1 = r'$\lambda_{1}=%.6f$' % (lam1, )
        eig_ax1.text(0.05, 0.9, textstr1, transform=eig_ax1.transAxes, fontsize=12, verticalalignment='top', linespacing=1.0)
        eig_ax1.set_xticks([])
        eig_ax1.set_yticks([])
        eig_ax1.set_title('Unstable Eigenvector')

        # second eigen value
        eig_ax2 = fig.add_subplot(gs_info[0,1])
        textstr2 = r'$\lambda_{2}=%.6f$' % (lam2, )
        eig_ax2.text(0.05, 0.9, textstr2, transform=eig_ax2.transAxes, fontsize=12, verticalalignment='top', linespacing=1.0)
        eig_ax2.set_xticks([])
        eig_ax2.set_yticks([])
        eig_ax2.set_title('Stable Eigenvector')

        # epsilon user entry
        eps_ax = fig.add_subplot(gs_info[1, 0]) 
        eps_box = TextBox(eps_ax, "Epsilon", textalignment="center", initial=epsilon, color='0.5', hovercolor="#4d4d4d")
        eps_box.on_submit(save_epsilon)

        # samples user entry
        samples_ax = fig.add_subplot(gs_info[1, 1]) 
        samples_box = TextBox(samples_ax, "Samples", textalignment="center", initial=samples, color='0.5', hovercolor="#4d4d4d")
        samples_box.on_submit(save_samples)
        
        # Tfactor user entry
        tfactor_ax = fig.add_subplot(gs_info[2, 0]) 
        tfactor_box = TextBox(tfactor_ax, "T_factor", textalignment="center", initial=T_factor, color='0.5', hovercolor="#4d4d4d")
        tfactor_box.on_submit(save_tfactor)
        
        # tsteps user entry
        tsteps_ax = fig.add_subplot(gs_info[2, 1]) 
        tsteps_box = TextBox(tsteps_ax, "tsteps", textalignment="center", initial=tsteps, color='0.5', hovercolor="#4d4d4d")
        tsteps_box.on_submit(save_tsteps)
        
        # recompute manifold button
        man_compute_ax = fig.add_subplot(gs_info[3,:] )
        man_compute_btn = Button(man_compute_ax, 'Recompute Manifolds', color="#2d2d2d" , hovercolor="#4d4d4d")
        man_compute_btn.label.set_color("white")
        man_compute_btn.on_clicked(refresh_manifolds)

        # manifold plot buttons
        unstable_btn_ax = fig.add_subplot(gs_info[4,0] )
        stable_btn_ax = fig.add_subplot(gs_info[4,1] )
        
        unstable_btn = Button(unstable_btn_ax, 'Unstable Manifold', color="#2d2d2d", hovercolor="#4d4d4d")
        unstable_btn.label.set_color("white")
        unstable_btn.on_clicked(plot_unstable_manifold)

        stable_btn = Button(stable_btn_ax, 'Stable Manifold', color="#2d2d2d" , hovercolor="#4d4d4d")
        stable_btn.label.set_color("white")
        stable_btn.on_clicked(plot_stable_manifold)
        
        
        #orbit_ax.axis('equal')
        
        # Reset the buttons
        info_buttons.clear()
        info_buttons.append([eps_box, samples_box, tfactor_box, tsteps_box, man_compute_btn, unstable_btn, stable_btn])


        # Clear old buttons and create new ones
        # buttons.clear()
        # for i, name in enumerate(["Orbit", "Phase Plots"]):
        #     button_ax = fig.add_axes([0.30 + i * 0.12, 0.95, 0.10, 0.03])
        #     btn = Button(button_ax, name, color="#2d2d2d", hovercolor="#4d4d4d")
        #     btn.label.set_color("white")
        #     if name == "Orbit":
        #         btn.on_clicked(lambda x: create_orbit_tab())
        #     elif name == "Phase Plots":
        #         btn.on_clicked(lambda x: create_phase_tab())
        #     #elif name == "Controls":
        #     #    btn.on_clicked(lambda x: create_controls_tab())
        #     #else:
        #     #    btn.on_clicked(lambda x: create_aero_tab())
        #     buttons.append(btn)
        
        plt.draw()
    

    def create_phase_tab():
        fig.clear()
        gs = fig.add_gridspec(1, 2, top=0.9)
        
        global buttons
        
        ## global figures

        
        buttons.clear()
        for i, name in enumerate(["Orbit", "Phase Plots"]):
            button_ax = fig.add_axes([0.30 + i * 0.12, 0.95, 0.10, 0.03])
            btn = Button(button_ax, name, color="#2d2d2d", hovercolor="#4d4d4d")
            btn.label.set_color("white")
            if name == "Orbit":
                btn.on_clicked(lambda x: create_orbit_tab())
            elif name == "Phase Plots":
                btn.on_clicked(lambda x: create_phase_tab())
            #elif name == "Controls":
            #    btn.on_clicked(lambda x: create_controls_tab())
            #else:
            #    btn.on_clicked(lambda x: create_aero_tab())
            buttons.append(btn)

        plt.draw()
        return
    create_orbit_tab()
        # add table or something with the two relevant eigenvectors
        # two buttons on bottom: Plot stable manifold and plot unstable manifold
        # will display the corresponding eigen values:
        # 

    return












# Create interactive dashboard
interactive_fig = plt.figure(figsize=(15, 10))
create_dashboard(interactive_fig)
plt.show()
