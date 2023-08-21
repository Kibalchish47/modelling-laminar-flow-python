# Modelling-Laminar-Flow-Python #
This project aims to solve the 2D Navier-Stokes equations using the finite difference method for single-phase laminar flow and verify results using the benchmark lid cavity test.

## 1. Introduction ##
Fluid flow can be understood using three main approaches: experimental, analytical, and numerical. 

The experimental approach involves conducting physical experiments in a laboratory to observe and measure fluid flow properties. 

The analytical approach involves deriving mathematical equations, such as the Navier-Stokes equations, to predict fluid behavior under simplified conditions. However, solving these equations analytically can be challenging. 

The numerical approach, known as Computational Fluid Dynamics (CFD), uses powerful computers to simulate and predict fluid flow by solving the equations numerically. CFD allows researchers to study a wide range of complex flow scenarios and optimize designs. It is particularly useful for preventing undesired fluid interactions, such as coffee spills, by analyzing flow dynamics and recommending modifications. 

This is the one that will be of our interest in this project, as I am neither capable nor willing to do experimental approach or to analytically solve the Navier-Stokes equations (if that's even possible).

## 2. Governing Equations ##
So, what is this set of equations that can completely describe how a fluid flows and where do they come from? Before answering the former question, let’s discuss the latter.

Consider a 2D box having a fixed volume in space. This is what we term the control volume.
![Alt text](<assets/Figure1 Control Volume.png>)

First, we will apply the principle of conservation of mass to the fluid in the control volume. For an incompressible fluid (most liquids), this means that whatever fluid enters the box must exit it. This is referred to as the equation of continuity in fluid mechanics.

Second, we will apply the principle of conservation of momentum to the control volume. This is slightly more abstract and complex compared to the previous case but eventually, this reduces to the incompressible Navier-Stokes equations.

## Note 1 ##
This project is mostly based on the project (with the same title as this one) proposed in this article: 
https://towardsdatascience.com/computational-fluid-dynamics-using-python-modeling-laminar-flow-272dad1ebec
