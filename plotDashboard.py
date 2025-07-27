import plot_tools
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
from matplotlib.widgets import Button, TextBox
#import cartopy.mpl.ticker as mticker
import warnings
from matplotlib import cm


import os
sol_dir = os.path.join(os.getcwd(),"OrbitSolutions")

## Extract state information from loaded solution set
orbit = 'Halo.csv' # Select orbit solution file
Sol = np.loadtxt(os.path.join(sol_dir, orbit), delimiter=',')
states = Sol[1:7,:]
tvec = Sol[0,:]
eps = 1*10**(-4)


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
        gs = fig.add_gridspec(4, 2, top=0.9)
        
        global buttons
        global info_buttons

        def save_epsilon(epsilon_value):
            global eps
            eps = epsilon_value
            print(eps)
            return 
            
        ## Normal 3D orbit plot
        ax = fig.add_subplot(gs[0:3, 0], projection='3d')
        plot_tools.Orbit3D(states, tvec, ax, args={'Frame':'Synodic'})

        
        ##eigen value selection stuff
        gs_info = gs[6].subgridspec(3, 2)

        eig_ax= fig.add_subplot(gs_info[0,0:2])
        lam1 = 9.5; lam2 = 0.2
        textstr = '\n'.join((
                    r'$\lambda_{1}=%.2f$ - Unstable Eigenvector' % (lam1, ),
                    r'$\lambda_{2}=%.2f$ - Stable Eigenvector' % (lam2, )))
        eig_ax.text(0.05, 0.9, textstr, transform=eig_ax.transAxes, fontsize=12, verticalalignment='top', linespacing=1.0)
        eig_ax.set_xticks([])
        eig_ax.set_yticks([])
        

        eps_box = TextBox(fig.add_subplot(gs_info[1, 0:2]), "Epsilon", textalignment="center", initial=0.0001, color='0.5')
        eps_box.on_submit(save_epsilon)


        btn1_ax = fig.add_subplot(gs_info[2,0] )
        btn2_ax = fig.add_subplot(gs_info[2,1] )
        btn_lam1 = Button(btn1_ax, 'Unstable Manifold', color="#2d2d2d", hovercolor="#4d4d4d")
        btn_lam1.label.set_color("white")
        btn_lam2 = Button(btn2_ax, 'Stable Manifold', color="#2d2d2d", hovercolor="#4d4d4d")
        btn_lam2.label.set_color("white")
        
        info_buttons.append([eps_box, btn_lam1, btn_lam2])

        # velocity_magnitude = np.sqrt(
        #     df["Velocity X (km/s)"] ** 2
        #     + df["Velocity Y (km/s)"] ** 2
        #     + df["Velocity Z (km/s)"] ** 2
        # )

        ## Create scatter plot with time-based coloring
        # velo_scatter = ax.scatter(
        #     velocity_magnitude,
        #     df["Altitude (km)"],
        #     c=df["Time (s)"],
        #     cmap="plasma",
        #     alpha=0.6,
        #     label="Velocity",
        # )

        # ax2=ax.twiny()
        # drag_plot = ax2.plot(
        #     df["Drag Force (N)"],
        #     df["Altitude (km)"],
        #     # c=df["Time (s)"],
        #     # cmap="plasma",
        #     # alpha=0.6,
        #     label="Drag Force",
        # )
        # # Add colorbar
        # cbar = plt.colorbar(velo_scatter)
        # cbar.set_label("Mission Elapsed Time [hours]")
        # cbar.ax.yaxis.set_major_formatter(
        #     plt.FuncFormatter(lambda x, p: f"{x / 3600:.1f}")
        # )

        # ax.set_xlabel("Velocity Magnitude [km/s]")
        # ax.set_ylabel("Altitude [km]")
        # ax2.set_xlabel("Drag Force [N]")
        # ax.set_title("Altitude vs Velocity (bot) and Drag ")
        # ax.grid(True)

        # lines, labels = ax.get_legend_handles_labels()
        # lines2, labels2 = ax2.get_legend_handles_labels()
        # ax2.legend(lines + lines2, labels + labels2, loc=0)
        # ax2.legend(loc="upper right")

        # Clear old buttons and create new ones
        buttons.clear()
        for i, name in enumerate(["Orbit", "Attitude", "Controls", "Aero"]):
            button_ax = fig.add_axes([0.30 + i * 0.12, 0.95, 0.10, 0.03])
            btn = Button(button_ax, name, color="#2d2d2d", hovercolor="#4d4d4d")
            btn.label.set_color("white")
            if name == "Orbit":
                btn.on_clicked(lambda x: create_orbit_tab())
            #elif name == "Attitude":
            #    btn.on_clicked(lambda x: create_attitude_tab())
            #elif name == "Controls":
            #    btn.on_clicked(lambda x: create_controls_tab())
            #else:
            #    btn.on_clicked(lambda x: create_aero_tab())
            buttons.append(btn)
        
        plt.draw()
        
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
