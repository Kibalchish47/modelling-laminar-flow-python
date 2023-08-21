# Modelling-Laminar-Flow-Python #
## The Project
This project aims to solve the 2D Navier-Stokes equations using the finite difference method for single-phase laminar flow and verify results using the benchmark lid cavity test (with NumPy vectorization for better performance).

The Navier-Stokes equations and the continuity equation are discretized using a second order finite difference method and solved numerically. To ensure a divergence free velocity at the end of each time step, the process of calculating velocity and pressure is as follows:
- Calculate intermediate starred velocities that may not necessarily be divergence free by solving the momentum equation without the pressure term. This is the predictor step.
- Differentitate the momentum equation (with the pressure term) and apply continuity to eliminate next time-step velocities. Thus, we are left with a Poisson equation for pressure in terms of the starred velocities. This step gives the pressure field for the next time-step.
- Calculate divergence-free velocities for the next time-step using the newly calculated pressure field and the starred velocities in a corrector step.

Table of Contents: 
1. [Introduction](#1-introduction)
2. [Governing Equations](#2-governing-equations)
3. [Numerical Methods](#3-numerical-methods)
4. [Code Organization](#4-code-organization)
5. [References](#references)

Let us now get a better feeling for the topic!

## 1. Introduction ##
Fluid flow can be understood using three main approaches: experimental, analytical, and numerical. 

The experimental approach involves conducting physical experiments in a laboratory to observe and measure fluid flow properties. 

The analytical approach involves deriving mathematical equations, such as the Navier-Stokes equations, to predict fluid behavior under simplified conditions. However, solving these equations analytically can be challenging. 

The numerical approach, known as Computational Fluid Dynamics (CFD), uses powerful computers to simulate and predict fluid flow by solving the equations numerically. CFD allows researchers to study a wide range of complex flow scenarios and optimize designs. It is particularly useful for preventing undesired fluid interactions, such as coffee spills, by analyzing flow dynamics and recommending modifications. 

This is the one that will be of our interest in this project, as I am neither capable nor willing to do experimental approach or to analytically solve the Navier-Stokes equations (if that's even possible).

## 2. Governing Equations ##
So, what is this set of equations that can completely describe how a fluid flows and where do they come from? Before answering the former question, let’s discuss the latter.

Consider a 2D box having a fixed volume in space. This is what we term the control volume.
![Alt text](<assets/Figure1 Control Volume.png> "Figure 1: Control Volume")

First, we will apply the principle of conservation of mass to the fluid in the control volume. For an incompressible fluid (most liquids), this means that whatever fluid enters the box must exit it. This is referred to as the equation of continuity in fluid mechanics.
$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$$

Second, we will apply the principle of conservation of momentum to the control volume. This is slightly more abstract and complex compared to the previous case but eventually, this reduces to the incompressible Navier-Stokes equations.
$$\frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial x} + \nu(\frac{\partial^{2} u}{\partial x^2} + \frac{\partial^{2} u}{\partial y^{2}})$$

$$\frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} = -\frac{1}{\rho}\frac{\partial p}{\partial y} + \nu(\frac{\partial^{2} v}{\partial x^2} + \frac{\partial^{2} v}{\partial y^{2}})$$

If we can solve these partial differential equations (PDEs) simultaneously after applying requisite boundary conditions, we will obtain the instantaneous velocities and pressure as a function of time, allowing us to predict how the fluid will flow. However, there is no analytical method to solve these equations (in their complete forms) without applying simplifying assumptions. Therefore, we resort to numerical techniques for solving these equations.

## 3. Numerical Methods ##
There exist a variety of different numerical methods for solving PDEs, each with its own set of caveats. The simplest method is the Finite Difference method wherein a low-order Taylor series approximation is used to convert the PDEs to a set of algebraic equations. An example is given below that shows how to convert first and second order derivatives to their finite difference approximations.

$$\frac{\partial u}{\partial t} \approx \frac{u(x+\Delta x) - u(x - \Delta x)}{2\Delta x} = \frac{\Delta u}{\Delta x}$$

$$\frac{\partial^{2} u}{\partial x^2} \approx \frac{u(x+\Delta x) - 2u(x) + u(x - \Delta x)}{2\Delta x} = \frac{\Delta^2 u}{\Delta x^2}$$

While this is not the best method to model fluid flow in all cases, we will proceed with it since it simplifies other aspects of modeling a crystallizer, which is the ultimate objective of this series of articles. For more rigorous numerical treatments, you may want to use the the Finite Volume or Finite Element methods.

## 4. Code Organization ##
The code is organized into three different files or scripts. 
- The first — “FlowPy.py” — contains the code for the solution of the PDEs using the finite difference method for a general set of inputs. 
- The inputs are provided to this script using the “FlowPy_Input.py” script which acts as a user interface. 
- Finally, the “FlowPy_Visualizer.py” script is used to animate the dynamics of the flow after running the simulation.

## References ##
This project is mostly based on the project (with the same title as this one) proposed in this article: 
https://towardsdatascience.com/computational-fluid-dynamics-using-python-modeling-laminar-flow-272dad1ebec

Other Important References:

- Barba, L. A., & Forsyth, G. F. (2018). CFD Python: the 12 steps to Navier-Stokes equations. Journal of Open Source Education, 2(16), 21.
- Owkes, M. (2020), A guide to writing your first CFD Solver
- Ghia, U. K. N. G., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. Journal of computational physics, 48(3), 387–411.