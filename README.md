# Modelling-Laminar-Flow-Python #
## The Project
This project aims to solve the 2D Navier-Stokes equations using the finite difference method for single-phase laminar flow and verify results using the benchmark lid cavity test (with NumPy vectorization for better performance).

The Navier-Stokes equations and the continuity equation are discretized using a second-order finite difference method and solved numerically. To ensure a divergence-free velocity at the end of each time step, the process of calculating velocity and pressure is as follows:
- Calculate intermediate starred velocities that may not necessarily be divergence-free by solving the momentum equation without the pressure term. This is the predictor step.
- Differentiate the momentum equation (with the pressure term) and apply continuity to eliminate next time-step velocities. Thus, we are left with a Poisson equation for pressure in terms of the starred velocities. This step gives the pressure field for the next time step.
- Calculate divergence-free velocities for the next time-step using the newly calculated pressure field and the starred velocities in a corrector step.

The Result:
![Alt text](Result/FluidFlowAnimation.gif)

Table of Contents: 
1. [Introduction](#1-introduction)
2. [Governing Equations](#2-governing-equations)
3. [Numerical Methods](#3-numerical-methods)
4. [Code Organization](#4-code-organization)
5. [Implementation](#5-implementation)
    - [Build Classes](#build-classes)
        - [Define (Mathematically) the Boundary](#define-mathematically-the-boundary)
        - [Define Domain enclosed by the Boundary](#define-domain-enclosed-by-the-boundary)
        - [Define Fluid Class](#define-fluid-class)
    - [Write Functions to Implement the Finite Difference Method](#write-functions-to-implement-the-finite-difference-method)
        - [Set boundary conditions for horizontal velocity](#set-boundary-conditions-for-horizontal-velocity)
        - [Set boundary conditions for vertical velocity](#set-boundary-conditions-for-vertical-velocity)
        - [Set boundary conditions for pressure](#set-boundary-conditions-for-pressure)
        - [Determine the time-step](#determine-the-time-step)
        - [Finite difference scheme](#finite-difference-scheme)
        - [Convenience Function](#convenience-function)
        - [I/O Functions](#io-functions)
    - [The Simulation User Interface: FlowPy_Input](#the-simulation-user-interface--flowpy_input)
        - [Imports](#flowpy_inputpy-imports)
        - [Define spatial, temporal, physical, and momentum parameters](#define-spacial-temporal-physical-and-momentum-parameters)
        - [Write the simulation loop](#write-the-simulation-loop)
    - [The Visualization Tool: FlowPy_Visualizer](#the-visualization-tool--flowpy_visualizer)
        - [Imports](#flowpy_visualizerpy-imports)
        - [Simulation inputs](#simulation-inputs)
        - [Text to array conversion](#text-to-array-conversion)
        - [Plot and save the animation](#plot-and-save-the-animation)
6. [Results](#6-results)
7. [Known issues](#7--known-issues)
8. [References](#references)

Let us now get a better feeling for the topic!

## 1. Introduction ##
Fluid flow can be understood using three main approaches: experimental, analytical, and numerical. 

The experimental approach involves conducting physical experiments in a laboratory to observe and measure fluid flow properties. 

The analytical approach involves deriving mathematical equations, such as the Navier-Stokes equations, to predict fluid behavior under simplified conditions. However, solving these equations analytically can be challenging. 

The numerical approach, known as Computational Fluid Dynamics (CFD), uses powerful computers to simulate and predict fluid flow by solving the equations numerically. CFD allows researchers to study a wide range of complex flow scenarios and optimize designs. It is particularly useful for preventing undesired fluid interactions, such as coffee spills, by analyzing flow dynamics and recommending modifications. 

This is the one that will be of interest in this project, as I am neither capable nor willing to do an experimental approach or to analytically solve the Navier-Stokes equations (if that's even possible).

## 2. Governing Equations ##
So, what is this set of equations that can completely describe how a fluid flows, and where do they come from? Before answering the former question, let’s discuss the latter.

Consider a 2D box having a fixed volume in space. This is what we term the control volume.
![Alt text](<Assets/Figure 1 Control Volume.png> "Figure 1: Control Volume")
Figure 1: Control Volume.

First, we will apply the principle of conservation of mass to the fluid in the control volume. For an incompressible fluid (most liquids), this means that whatever fluid enters the box must exit it. This is referred to as the equation of continuity in fluid mechanics.

$$
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0
$$

Second, we will apply the principle of conservation of momentum to the control volume. This is slightly more abstract and complex compared to the previous case but eventually, this reduces to the incompressible Navier-Stokes equations.

$$
\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu(\frac{\partial^{2} u}{\partial x^2} + \frac{\partial^{2} u}{\partial y^{2}})
$$

$$
\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial y} + \nu(\frac{\partial^{2} v}{\partial x^2} + \frac{\partial^{2} v}{\partial y^{2}})
$$

If we can solve these partial differential equations (PDEs) simultaneously after applying requisite boundary conditions, we will obtain the instantaneous velocities and pressure as a function of time, allowing us to predict how the fluid will flow. However, there is no analytical method to solve these equations (in their complete forms) without applying simplifying assumptions (that we know of). Therefore, we resort to numerical techniques for solving these equations. If you find one, contact the Clay Mathematics Institute to claim your 1,000,000$.

## 3. Numerical Methods ##
There exist a variety of different numerical methods for solving PDEs, each with its own set of caveats. The simplest method is the Finite Difference method wherein a low-order Taylor series approximation is used to convert the PDEs to a set of algebraic equations. An example is given below that shows how to convert first and second-order derivatives to their finite difference approximations.

$$
\frac{\partial u}{\partial t} \approx \frac{u(x+\Delta x) - u(x - \Delta x)}{2\Delta x} = \frac{\Delta u}{\Delta x}
$$

$$
\frac{\partial^{2} u}{\partial x^2} \approx \frac{u(x+\Delta x) - 2u(x) + u(x - \Delta x)}{2\Delta x} = \frac{\Delta^2 u}{\Delta x^2}
$$

While this is not the best method to model fluid flow in all cases, we will proceed with it since it simplifies other aspects of modeling a crystallizer. For more rigorous numerical treatments, you may want to use the Finite Volume or Finite Element methods (maybe I'll make a separate project utilizing those methods).

## 4. Code Organization ##
The code is organized into three different files or scripts. 
- The first — “FlowPy.py” — contains the code for the solution of the PDEs using the finite difference method for a general set of inputs. 
- The inputs are provided to this script using the “FlowPy_Input.py” script which acts as a user interface. 
- Finally, the “FlowPy_Visualizer.py” script is used to animate the dynamics of the flow after running the simulation.

## 5. Implementation
### Build Classes
We start by making classes for specific properties of the problem, starting with the boundary condition. The PDEs are solved by applying certain boundary conditions which indicate how the fluid will behave at the boundaries of the domain. For example, fluid flowing through a pipe will have a wall with zero fluid velocity and an entry as well as an exit with some specified flow velocity.
#### Define (Mathematically) the Boundary 
Mathematically, boundary conditions can be expressed in two forms — Dirichlet and Neumann boundaries. The former specifies a value of the dependent variable at the boundary whereas the latter specifies a value for the derivative of the dependent variable at the boundary. Therefore, we make a Boundary class that has two properties — type and value.
See FlowPy.py (lines 7 to 19).
#### Define Domain enclosed by the Boundary
Next, the domain enclosed by the boundary (like the inside of a pipe) is represented using a 2D mesh or grid, and the values of dependent variables are calculated at the center of boxes in the grid (for pressure) or at the faces of the boxes (for velocities). This is referred to as a staggered grid approach. To represent the mesh, we create a class called Space. The method CreateMesh creates a matrix of given size for the dependent variables and the SetDeltas method calculates the values of the differential lengths based on the specified length and breadth of the domain.   
See FlowPy.py (lines 19 to 71).

#### Define Fluid Class
Lastly, we create a class Fluid to represent the properties of the fluid—like density (rho) and viscosity (mu).
See FlowPy.py (lines 71 to 80).

### Write Functions to Implement the Finite Difference Method
As in the previous section, we start by writing functions to implement boundary conditions for the horizontal velocity (u), vertical velocity (v) and pressure (p) at the left, right, top and bottom boundaries of the 2D domain. This function will accept the objects of the Space and Boundary classes and set boundary conditions according to the attributes of those objects. For example, if a Boundary object with type Dirichlet and value 0 is passed as the left boundary object, the function will set that condition on the left boundary.
#### Set boundary conditions for horizontal velocity
See FlowPy.py (lines 89 to 108).
#### Set boundary conditions for vertical velocity
See FlowPy.py (lines 71 to 80).
#### Set boundary conditions for pressure
See FlowPy.py (lines 71 to 80).
#### Determine the time-step
Before we write the finite difference functions, we need to determine a time step to advance the simulation by. To ensure the convergence of finite difference methods, an upper bound on the time-step is provided by the Courant–Friedrichs–Lewy (CFL) criterion which is set as the time-step for the simulation using the SetTimeStep function. Adhering to the CFL criterion ensures that information propagated in a time step is not farther than the distance between two mesh elements.
See FlowPy.py (lines 153 to 168).
#### Finite difference scheme
Having determined the time step, we are now ready to implement the finite difference scheme. To solve the equation of continuity and the Navier-Stokes equations simultaneously, we use a predictor-corrector scheme involving the following steps (for more information refer to this guide): https://www.montana.edu/mowkes/research/source-codes/GuideToCFD_2020_02_28_v2.pdf

- Calculate starred velocities (u* and v*) from initial velocities without the effect of pressure.

$$ 
u^{*}(t) = u(t) + \Delta t \left[-u(t)\frac{\Delta u(t)}{\Delta x} - v(t)\frac{\Delta u(t)}{\Delta y} 
+ \nu \left (\frac{\Delta^{2} u(t)}{\Delta x^{2}} + \frac{\Delta^{2} u(t)}{\Delta y^{2}} \right) \right]
$$

- Iteratively solve the pressure Poisson equation using the starred velocities.

$$ 
\frac{\Delta^{2} p(t + \Delta t)}{\Delta x^{2}} + \frac{\Delta^{2} p(t + \Delta t)}{\Delta y^{2}} = -\frac{\rho}{\Delta t} \left(\frac{\Delta u^{*}(t)}{\Delta x} + \frac{\Delta u^{*}(t)}{\Delta y} \right)
$$

- Calculate the velocities for the next time-step from the pressure and starred velocities.

$$ 
u(t + \Delta t) = u^{*}(t) + \Delta t \left(-\frac{1}{\rho} \frac{\Delta p}{\Delta x} \right)
$$

We define three different functions to carry out each of these three steps.
See FlowPy.py (lines 169 to 288). 
#### Convenience Function
Further, a convenience function is defined to save the velocities and pressures inside the boundaries to new variables, which can then be written to text files.
See FlowPy.py (lines 288 to 295).
#### I/O Functions
Finally, we define two functions for I/O purposes — MakeResultDirectory to make a directory called “Result” to store the text files and WriteToFile to save the values of the variables to a text file every few iterations (specified using the interval argument).
See FlowPy.py (lines 295 to 325).

### The Simulation User Interface — FlowPy_Input
This section is shorter than the previous one — most of the heavy lifting has been done, we just need to make use of all the defined classes and functions to run the simulation now!

As an example, inputs relevant to the Lid Cavity Test (at Reynolds Number=400) are entered in this tutorial. In this test, fluid is kept in a 2D box with three rigid walls and the fourth wall (or the lid) is moved at a uniform velocity. Once a steady state is reached, statistics of the developed flow field can be compared to a benchmark.

![Alt text](<Assets/Figure 2 Lid Cavity Problem set-up.png> "Figure 2: Lid Cavity Problem set-up")
Figure 2: Lid Cavity Problem set-up
#### FlowPy_Input.py Imports
First, we import the required modules and this now includes all the things that we have defined in FlowPy.py
See FlowPy_Input.py (lines 1 to 6). 
#### Define spatial, temporal, physical, and momentum parameters
We begin by specifying input variables describing the domain and then creating a Space object with these variables.
Next, the density and viscosity of the fluid are specified, and an object of the class Fluid is created.
Third, we create Boundary objects to set velocity and pressure boundary conditions.
Finally, simulation parameters and flags are specified that control simulation time, saving text files, and so forth.
See FlowPy_Input.py (lines 13 to 48). 

#### Write the simulation loop
Now, we can write the loop to run the simulation. The general procedure is as follows. Until the simulation time is completed, do the following in every iteration:
- Set the time step according to the CFL number criterion
- Set boundary conditions
- Calculate starred velocities
- Solve the pressure Poisson equation to get the pressure field
- Determine velocities at the next time-step
- Write results to file (if file flag is 1)
- Advance time by a value equal to the time step
See FlowPy_Input.py (lines 48 to 101)
Having reached here, we are now ready to run the simulation for any generalized set of inputs. There’s just one piece of the puzzle left — a visualization tool.

### The Visualization Tool — FlowPy_Visualizer
The text files that are generated after running the simulation contain raw numbers that may not provide a physical picture of the fluid flow by themselves. However, a simple, animated contour plot can be used to combine the three variables — horizontal velocity, vertical velocity, and pressure — and show their time evolution in an intuitive manner.

#### FlowPy_Visualizer.py Imports
As before, first, import the required modules. Particularly, we will require the matplotlib.animation module to record the animation.
See FlowPy_Visualizer.py (lines 1 to 7). 
#### Simulation inputs
To ensure that arrays of appropriate sizes are created, simulation inputs pertaining to the computational domain need to be entered.
See FlowPy_Visualizer.py (lines 13 to 21). 
#### Text to array conversion
Before moving to plotting, the text files that were saved during the simulation have to be imported as arrays. To do so, we first go through the Result directory, store all the filenames, and determine the total number of files as well as the printing interval.
Next, we define a function that can import a text file — based on a provided iteration — into an array using the loadtxt function in numpy.
See FlowPy_Visualizer.py (lines 21 to 75).
#### Plot and save the animation
It’s time to start making the plot! Before animating the figure, it’s a good idea to make an initial plot (for the zeroth iteration) so that the figure dimensions, axes, color bar and so on can be fixed. Also, it’s a good idea to make the stream plot with fewer grid points (here, 10) to make the arrows distinguishable.
To animate this plot further, the FuncAnimation function from matplotlib.animation will come in handy. All it needs is a function that can create a plot for a supplied value of the iteration. We define such a function called animate.
Finally, save the animation and watch some fluids dance around on your computer!
See FlowPy_Visualizer.py (lines 75 to 149).

## 6. Results ##
The animated contour and stream plot for the lid cavity benchmark at Re=400 is shown below. It shows the formation of a vortex at the center as the simulation progresses and ultimately, a transition to a steady state.

![Alt text](Result/FluidFlowAnimation.gif, "Resulting plot")

A quantitative comparison of the statistics of the steady flow with the results of Ghia et al. (1982) is also performed. Specifically, the horizontal velocities along a vertical line passing through the center of the cavity and vice versa are compared with the simulation results from the paper. The results show reasonable agreement. Deviations can be attributed to the lower accuracy of the finite difference scheme and smaller grid size.

![Alt text](Assets/Figure%203.webp, "Figure 3: Benchmark 1. The blue line represents simulation results and the red points represent results from Ghia et al. (1982).")
Figure 3: Benchmark 1. The blue line represents simulation results and the red points represent results from Ghia et al. (1982).

![Alt text](Assets/Figure%204.webp, "Figure 4: Benchmark 2. The blue line represents simulation results and the red points represent results from Ghia et al. (1982).")
Figure 4: Benchmark 2. The blue line represents simulation results and the red points represent results from Ghia et al. (1982). 

While this tutorial only includes the simulation of the lid cavity test, you can try playing around with the inputs and boundary conditions to model a variety of different single-phase flow problems, like Poiseuille flow in a pipe.

With the creation and validation of FlowPy, we can move to the next step in the modeling of a crystallizer — the addition of heat and mass transfer to the solver, which will be covered in the next article.
## 7. Known Issues ##
There aren't many issues with the project, but here are some I encountered (and you need to be aware of too): 
- Some parts of the code are extremely bulky, especially the ones involving NumPy vectorization (and I have no clue how to simplify it)
- It took a while to make the simulation save properly. If you have any issues with saving the simulation, check out this tutorial: https://holypython.com/how-to-save-matplotlib-animations-the-ultimate-guide/?expand_article=1

## 8. References ##
- This project is purely based on the project (with the same title as this one) proposed in this article: https://towardsdatascience.com/computational-fluid-dynamics-using-python-modeling-laminar-flow-272dad1ebec. The article's code is also available on GitHub: https://github.com/gauravsdeshmukh/FlowPy/tree/master.
- https://www.montana.edu/mowkes/research/source-codes/GuideToCFD_2020_02_28_v2.pdf
- https://www.cfd-online.com/Wiki/Lid-driven_cavity_problem

Other Important References:

- Barba, L. A., & Forsyth, G. F. (2018). CFD Python: the 12 steps to Navier-Stokes equations. Journal of Open Source Education, 2(16), 21.
- Owkes, M. (2020), A guide to writing your first CFD Solver
- Ghia, U. K. N. G., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. Journal of computational physics, 48(3), 387–411. 
