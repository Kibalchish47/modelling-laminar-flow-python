# Basic imports (yeah it's not much)
import numpy as np
import numba as nb
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
# ---------------------------------------------- 
# Mathematically, boundary conditions can be expressed in two forms — Dirichlet and Neumann boundaries. 
# The former specifies a value of the dependent variable at the boundary 
# whereas the latter specifies a value for the derivative of the dependent variable at the boundary.
# Therefore, we make a Boundary class that has two properties — type and value.
class Boundary:
    def __init__(self, boundary_type, boundary_value):
        self.DefineBoundary(boundary_type, boundary_value)
        
    def DefineBoundary(self, boundary_type, boundary_value):
        self.type = boundary_type
        self.value = boundary_value
# ---------------------------------------------- 
# Next, the domain enclosed by the boundary (like the inside of a pipe) is represented using a 2D mesh or grid 
# and the values of dependent variables are calculated at the center of boxes in the grid (for pressure) 
# or at the faces of the boxes (for velocities). This is referred to as a staggered grid approach. 
# To represent the mesh, we create a class called Space. The method CreateMesh creates a matrix of given size 
# for the dependent variables and the SetDeltas method calculates the values of the differential lengths based 
# on the specified length and breadth of the domain.   
    
class Space:
    def __init__(self):
        pass
    def CreateMesh(self, rowpts, colpts):
        #Domain gridpoints
        self.rowpts = rowpts
        self.colpts = colpts
        
        #Velocity matrices
        self.u = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.v = np.zeros((self.rowpts + 2, self.colpts + 2))
        
        self.u_star = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.v_star = np.zeros((self.rowpts + 2, self.colpts + 2))
        
        self.u_next = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.v_next = np.zeros((self.rowpts + 2, self.colpts + 2))
        
        self.u_c = np.zeros((self.rowpts, self.colpts))
        self.v_c = np.zeros((self.rowpts, self.colpts))
        
        #Pressure matrices
        self.p = np.zeros((self.rowpts + 2, self.colpts + 2))
        self.p_c = np.zeros((self.rowpts, self.colpts))

        #Set default source term
        self.SetSourceTerm()        
        
    def SetDeltas(self,breadth, length):
        self.dx = length / (self.colpts - 1)
        self.dy = breadth / (self.rowpts - 1)
        
    def SetInitialU(self, U):
        self.u = U * self.u
        
    def SetInitialV(self, V):
        self.v = V * self.v
        
    def SetInitialP(self, P):
        self.p = P * self.p
        
    def SetSourceTerm(self, S_x = 0, S_y = 0):
        self.S_x = S_x
        self.S_y = S_y
# ---------------------------------------------- 
# Lastly, we create a class Fluid to represent the properties of the fluid — like density (rho) and viscosity (mu).
class Fluid:
    def __init__(self, rho, mu):
        self.SetFluidProperties(rho, mu)
    
    def SetFluidProperties(self, rho, mu):
        self.rho = rho
        self.mu = mu
# ---------------------------------------------- 
# As in the previous section, we start by writing functions to implement boundary conditions 
# for the horizontal velocity (u), vertical velocity (v) and pressure (p) at the left, right, top and bottom boundaries 
# of the 2D domain. This function will accept the objects of the Space and Boundary classes and set boundary conditions 
# according to the attributes of those objects. For example, if a Boundary object with type Dirichlet and value 0 
# is passed as the left boundary object, the function will set that condition on the left boundary.
# Note: The arguments to the function are all objects of our defined classes
# ---------------------------------------------- 
#### Set boundary conditions for horizontal velocity
def SetUBoundary(space, left, right, top, bottom):
    if(left.type == "D"):
        space.u[:, 0] = left.value
    elif(left.type == "N"):
        space.u[:, 0] = -left.value * space.dx + space.u[:, 1]
    
    if(right.type == "D"):
        space.u[:, -1] = right.value
    elif(right.type == "N"):
        space.u[:, -1] = right.value * space.dx + space.u[:, -2]
        
    if(top.type == "D"):
        space.u[-1, :] = 2 * top.value - space.u[-2,:]
    elif(top.type=="N"):
        space.u[-1, :] = -top.value * space.dy + space.u[-2, :]
     
    if(bottom.type == "D"):
        space.u[0, :] = 2 * bottom.value - space.u[1,:]
    elif(bottom.type == "N"):
        space.u[0, :] = bottom.value * space.dy + space.u[1,:]
# ----------------------------------------------      
#### Set boundary conditions for vertical velocity
def SetVBoundary(space, left, right, top, bottom):
    if(left.type == "D"):
        space.v[:, 0] = 2 * left.value - space.v[:, 1]
    elif(left.type == "N"):
        space.v[:, 0] = -left.value * space.dx + space.v[:, 1]
    
    if(right.type == "D"):
        space.v[:, -1] = 2 * right.value - space.v[:, -2]
    elif(right.type == "N"):
        space.v[:, -1] = right.value * space.dx + space.v[:, -2]
        
    if(top.type == "D"):
        space.v[-1, :] = top.value
    elif(top.type == "N"):
        space.v[-1, :] = -top.value * space.dy + space.v[-2, :]
     
    if(bottom.type == "D"):
        space.v[0, :] = bottom.value
    elif(bottom.type == "N"):
        space.v[0, :] = bottom.value * space.dy + space.v[1, :]
# ---------------------------------------------- 
#### Set boundary conditions for pressure
def SetPBoundary(space, left, right, top, bottom):
    if(left.type == "D"):
        space.p[:, 0] = left.value
    elif(left.type == "N"):
        space.p[:, 0] = -left.value * space.dx + space.p[:, 1]
    
    if(right.type == "D"):
        space.p[1, -1] = right.value
    elif(right.type == "N"):
        space.p[:, -1] = right.value * space.dx + space.p[:, -2]
        
    if(top.type == "D"):
        space.p[-1, :] = top.value
    elif(top.type == "N"):
        space.p[-1, :] = -top.value * space.dy + space.p[-2, :]
     
    if(bottom.type == "D"):
        space.p[0, :] = bottom.value
    elif(bottom.type == "N"):
        space.p[0, :] = bottom.value * space.dy + space.p[1, :]
# ---------------------------------------------- 
## To ensure the convergence of finite difference methods, an upper bound on the time-step is 
## provided by the Courant–Friedrichs–Lewy (CFL) criterion which is set as the time-step for 
## the simulation using the SetTimeStep function. 
## Adhering to the CFL criterion ensures that information propagated in a time-step is 
## not farther than the distance between two mesh elements.

def SetTimeStep(CFL, space, fluid):
    with np.errstate(divide = 'ignore'):
        dt=CFL/np.sum([np.amax(space.u) / space.dx,\
                        np.amax(space.v)/ space.dy])
        
    # Escape condition if dt is infinity due to zero velocity initially
    if np.isinf(dt):
        dt = CFL * (space.dx + space.dy)
    space.dt = dt
# ---------------------------------------------- 
# Having determined the time-step, we are now ready to implement the finite difference scheme.
#### We define three different functions to carry out each of these three steps.
#### The first function is used to get starred velocities from u and v at timestep t without the effect of pressure
def GetStarredVelocities(space, fluid):
    # Save object attributes as local variable (with explicit typing for improved readability)
    # from the space and fluid objects, such as the number of rows 
    # and columns in the computational domain (rows and cols), the velocity components (u and v), 
    # and various parameters like the grid spacing (dx and dy), time step size (dt), fluid density (rho), and fluid viscosity (mu).
    rows = int(space.rowpts)
    cols = int(space.colpts)
    u = space.u.astype(float, copy=False)
    v = space.v.astype(float, copy=False)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    S_x = float(space.S_x)
    S_y = float(space.S_y)
    rho = float(fluid.rho)
    mu = float(fluid.mu)
    
    #Copy u and v to new variables u_star and v_star
    u_star = u.copy()
    v_star = v.copy()
    
    # Calculate derivatives of u and v using the finite difference scheme. 
    # Numpy vectorization saves us from using slower for loops to go over each element in the u and v matrices
    u1_y = (u[2:, 1:cols+1] - u[0:rows, 1:cols+1]) / (2*dy)
    u1_x = (u[1:rows+1, 2:] - u[1:rows+1, 0:cols]) / (2*dx)
    u2_y = (u[2:, 1:cols+1] - 2*u[1:rows+1, 1:cols+1] + u[0:rows, 1:cols+1]) / (dy**2)
    u2_x = (u[1:rows+1, 2:] - 2*u[1:rows+1, 1:cols+1] + u[1:rows+1, 0:cols]) / (dx**2)
    v_face = (v[1:rows+1, 1:cols+1] + v[1:rows+1, 0:cols] + v[2:, 1:cols+1] + v[2:, 0:cols]) / 4
    u_star[1:rows+1, 1:cols+1] = u[1:rows+1, 1:cols+1] - dt*(u[1:rows+1, 1:cols+1]*u1_x + v_face*u1_y) \
        + (dt*(mu/rho)*(u2_x+u2_y)) + (dt*S_x)   

    v1_y = (v[2:, 1:cols+1] - v[0:rows, 1:cols+1]) / (2*dy)
    v1_x = (v[1:rows+1, 2:] - v[1:rows+1, 0:cols]) / (2*dx)
    v2_y = (v[2:, 1:cols+1] - 2*v[1:rows+1, 1:cols+1] + v[0:rows, 1:cols+1]) / (dy**2)
    v2_x = (v[1:rows+1, 2:] - 2*v[1:rows+1, 1:cols+1] + v[1:rows+1, 0:cols]) / (dx**2)
    u_face = (u[1:rows+1, 1:cols+1] + u[1:rows+1, 2:] + u[0:rows, 1:cols+1] + u[0:rows, 2:]) / 4
    v_star[1:rows+1, 1:cols+1] = v[1:rows+1, 1:cols+1] - dt*(u_face*v1_x+v[1:rows+1, 1:cols+1]*v1_y) \
        +(dt*(mu/rho)*(v2_x + v2_y))+(dt*S_y)
    
    #Save the calculated starred velocities to the space object 
    space.u_star = u_star.copy()
    space.v_star = v_star.copy()    

#### The second function is used to iteratively solve the pressure Possion equation from the starred velocities 
#### to calculate pressure at t+delta_t
def SolvePressurePoisson(space, fluid, left, right, top, bottom):
    # Save object attributes as local variable with explicitly typing for improved readability
    rows = int(space.rowpts)
    cols = int(space.colpts)
    u_star = space.u_star.astype(float, copy=False)
    v_star = space.v_star.astype(float, copy=False)
    p = space.p.astype(float, copy=False)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    rho = float(fluid.rho)
    factor = 1 / (2/dx**2 + 2/dy**2)
    
    # Define initial error and tolerance for convergence (error > tol necessary initially)
    error = 1
    tol = 1e-3

    # Evaluate derivative of starred velocities
    ustar1_x = (u_star[1:rows+1, 2:] - u_star[1:rows+1, 0:cols]) / (2*dx)
    vstar1_y = (v_star[2:, 1:cols+1] - v_star[0:rows, 1:cols+1]) / (2*dy)

    # Continue iterative solution until error becomes smaller than tolerance
    i = 0
    while(error > tol):
        i += 1
        
        #Save current pressure as p_old
        p_old=p.astype(float, copy=True)
        
        #Evaluate second derivative of pressure from p_old
        p2_xy = (p_old[2:,1:cols+1] + p_old[0:rows,1:cols+1]) / dy**2 + (p_old[1:rows+1,2:] + p_old[1:rows+1,0:cols]) / dx**2
        
        #Calculate new pressure 
        p[1:rows+1,1:cols+1] = (p2_xy) * factor - (rho*factor/dt) * (ustar1_x+vstar1_y)
        
        #Find maximum error between old and new pressure matrices
        error = np.amax(abs(p-p_old))
        
        #Apply pressure boundary conditions
        SetPBoundary(space, left, right, top, bottom)
        
        #Escape condition in case solution does not converge after 500 iterations
        if(i > 500):
            tol *= 10
# ---------------------------------------------- 
#### The third function is used to calculate the velocities at timestep t+delta_t 
#### using the pressure at t+delta_t and starred velocities
def SolveMomentumEquation(space, fluid):
    #Save object attributes as local variable with explicitly typing for improved readability
    rows = int(space.rowpts)
    cols = int(space.colpts)
    u_star = space.u_star.astype(float, copy=False)
    v_star = space.v_star.astype(float, copy=False)
    p = space.p.astype(float, copy=False)
    dx = float(space.dx)
    dy = float(space.dy)
    dt = float(space.dt)
    rho = float(fluid.rho)
    u = space.u.astype(float, copy=False)
    v = space.v.astype(float, copy=False)

    #Evaluate first derivative of pressure in x direction
    p1_x = (p[1:rows+1, 2:] - p[1:rows+1, 0:cols]) / (2*dx)
    #Calculate u at next timestep
    u[1:rows+1, 1:cols+1] = u_star[1:rows+1, 1:cols+1] - (dt/rho) * p1_x

    #Evaluate first derivative of pressure in y direction
    p1_y = (p[2:, 1:cols+1] - p[0:rows, 1:cols+1]) / (2*dy)
    #Calculate v at next timestep
    v[1:rows+1, 1:cols+1] = v_star[1:rows+1, 1:cols+1] - (dt/rho) * p1_y
# ---------------------------------------------- 
#### Convenience function to save the velocities and pressures inside the boundaries to new variables,
#### which can then be written to text files.
def SetCentrePUV(space):
    space.p_c = space.p[1:-1, 1:-1]
    space.u_c = space.u[1:-1, 1:-1]
    space.v_c = space.v[1:-1, 1:-1]
# ---------------------------------------------- 
# Finally, we define two functions for I/O purposes: 
#### - MakeResultDirectory to make a directory called “Result to store the text files
#### - WriteToFile to save the values of the variables to a text file every few iterations 
####   (specified using the interval argument).
def MakeResultDirectory(wipe=False):
    #Get path to the Result directory
    cwdir = os.getcwd()
    dir_path = os.path.join(cwdir, "Result")
    #If directory does not exist, make it
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path, exist_ok=True)
    else:
        #If wipe is True, remove files present in the directory
        if wipe:
            os.chdir(dir_path)
            filelist = os.listdir()
            for file in filelist:
                os.remove(file)
    os.chdir(cwdir)           
    
def WriteToFile(space, iteration, interval):
    if(iteration % interval == 0):
        dir_path = os.path.join(os.getcwd(), "Result")
        filename = "PUV{0}.txt".format(iteration)
        path = os.path.join(dir_path, filename)
        with open(path, "w") as f:
            for i in range(space.rowpts):
                for j in range(space.colpts):
                    f.write("{}\t{}\t{}\n".format(space.p_c[i,j],space.u_c[i,j],space.v_c[i,j]))
# ---------------------------------------------- 
