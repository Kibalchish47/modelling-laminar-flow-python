# Imports
import numpy as np
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
mpl.rcParams['animation.ffmpeg_path'] = r'C:\Users\kibal\OneDrive\Документы'
# ----------------------------------------------
# The text files that are generated after running the simulation contain raw numbers 
# that may not provide a physical picture of the fluid flow by themselves. 
# However, a simple, animated contour plot can be used to combine the three variables — horizontal velocity, 
# vertical velocity and pressure — and show their time evolution in an intuitive manner.
# ----------------------------------------------
# To ensure that arrays of appropriate sizes are created, 
# simulation inputs pertaining to the computational domain need to be entered.
#### Simulation inputs
rowpts = 257
colpts = 257
length = 4
breadth = 4
# ----------------------------------------------
# Before moving to plotting, the text files that were saved during the simulation have 
# to be imported as arrays. To do so, we first go through the Result directory, store all the filenames, 
# and determine the total number of files as well as the printing interval.
#### Go to the Result directory
cwdir = os.getcwd()
dir_path = os.path.join(cwdir,"Result")
os.chdir(dir_path)

# Go through files in the directory and store filenames
filenames = []
iterations = []
for root,dirs,files in os.walk(dir_path):
    for datafile in files:
        if "PUV" in datafile:
            filenames.append(datafile)
            no_ext_file=datafile.replace(".txt", "").strip()
            iter_no=int(no_ext_file.split("V")[-1])
            iterations.append(iter_no)

# Discern the final iteration and interval
initial_iter = np.amin(iterations)            
final_iter = np.amax(iterations)
inter = (final_iter - initial_iter) / (len(iterations)-1)
number_of_frames = len(iterations)
sorted_iterations = np.sort(iterations)
# ----------------------------------------------
# Next, we define a function that can import a text file — based on a provided iteration 
# — into an array using the loadtxt function in numpy.

def read_datafile(iteration):
    #Set filename and path according to given iteration
    filename = "PUV{0}.txt".format(iteration)
    filepath = os.path.join(dir_path, filename)
    
    #Load text file as numpy array
    arr = np.loadtxt(filepath, delimiter="\t")
    rows, cols = arr.shape
    
    #Define empty arrays for pressure and velocities
    p_p = np.zeros((rowpts, colpts))
    u_p = np.zeros((rowpts, colpts))
    v_p = np.zeros((rowpts, colpts))
    
    #Organize imported array into variables
    p_arr = arr[:,0]
    u_arr = arr[:,1]
    v_arr = arr[:,2]
    
    #Reshape 1D data into 2D
    p_p = p_arr.reshape((rowpts,colpts))
    u_p = u_arr.reshape((rowpts,colpts))
    v_p = v_arr.reshape((rowpts,colpts))
    
    return p_p, u_p, v_p
# ---------------------------------------------- 
# It’s time to start making the plot! Before animating the figure, it’s a good idea 
# to make an initial plot (for the zeroth iteration) so that the figure dimensions, axes, 
# color bar and so on can be fixed. Also, it’s a good idea to make the stream plot with 
# fewer grid points (in this article, 10) to make the arrows distinguishable.

# Create mesh for X and Y inputs to the figure
x = np.linspace(0,length,colpts)
y = np.linspace(0,breadth,rowpts)
[X,Y] = np.meshgrid(x,y)

# Determine indexing for stream plot (10 points only)
index_cut_x = int(colpts/10)
index_cut_y = int(rowpts/10)

# Create blank figure
fig = plt.figure(figsize=(16, 8))
ax = plt.axes(xlim=(0,length), ylim=(0, breadth))

# Create initial contour and stream plot as well as color bar
p_p, u_p, v_p = read_datafile(0)
ax.set_xlim([0, length])
ax.set_ylim([0, breadth])
ax.set_xlabel("$x$", fontsize=12)
ax.set_ylabel("$y$", fontsize=12)
ax.set_title("Frame No: 0")
cont = ax.contourf(X,Y, p_p)
stream = ax.streamplot(X[::index_cut_y,::index_cut_x],\
                       Y[::index_cut_y,::index_cut_x],\
                       u_p[::index_cut_y,::index_cut_x],\
                       v_p[::index_cut_y,::index_cut_x],\
                       color="k")
fig.colorbar(cont)
fig.tight_layout()
# ---------------------------------------------- 
# To animate this plot further, the FuncAnimation function from matplotlib.animation will come in handy. 
# All it needs is a function that can create a plot for a supplied value of the iteration. 
# We define such a function called animate.
def animate(i):
    #Print frames left to be added to the animation
    sys.stdout.write("\rFrames remaining: {0:03d}".format(len(sorted_iterations)-i))
    sys.stdout.flush()
    
    #Get iterations in a sequential manner through sorted_iterations
    iteration=sorted_iterations[i]
    
    #Use the read_datafile function to get pressure and velocities
    p_p,u_p,v_p=read_datafile(iteration)
    
    #Clear previous plot and make contour and stream plots for current iteration
    ax.clear()
    ax.set_xlim([0, length])
    ax.set_ylim([0, breadth])
    ax.set_xlabel("$x$", fontsize=12)
    ax.set_ylabel("$y$", fontsize=12)
    ax.set_title("Frame No: {0}".format(i))
    cont = ax.contourf(X, Y, p_p)
    stream = ax.streamplot(X[::index_cut_y,::index_cut_x],\
                            Y[::index_cut_y,::index_cut_x],\
                            u_p[::index_cut_y,::index_cut_x],\
                            v_p[::index_cut_y,::index_cut_x],\
                            color="k")
    return cont, stream
# ---------------------------------------------- 
# Finally, it’s time to save the animation and watch some fluids dance around on your computer!
print("######## Making FlowPy Animation ########")
print("#########################################")
anim = animation.FuncAnimation(fig, animate, frames=number_of_frames, interval=50, blit=False)
movie_path = os.path.join(dir_path,"FluidFlowAnimation.mp4")
writergif = animation.PillowWriter(fps=30) 
anim.save("FluidFlowAnimation.gif")
print("\nAnimation saved as FluidFlowAnimation.gif in Result")
# ---------------------------------------------- 