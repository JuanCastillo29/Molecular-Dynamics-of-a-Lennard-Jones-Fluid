# Molecular-Dynamics-of-a-Lennard-Jones-Fluid
This repository contains a basic code in C++ for simulate a Lennard-Jones fluid developed for the Molecular Modelling course of the MSc 'Physics of Complex Systems and biophysics' at Barcelona Univerity. The project was divided in two different exercises so I divided the code into two directories.

## Whitout Thermal Reservoir

In this section, I performed the simulation in the (E, V, N) ensemble so no thermal bath is implemented in the main section of the code. This exercise allows to compare the results obtained with different integration algorithms (mainly Euler and Velocity Verlet).

The production code is called 'Codigo.cpp'. It takes a value of the density and time step (in reduced units) and produces four files:
  - 'Kinetics.txt'. Five column file with time (in time steps), momentum in x, momentum in y, momentum in z and total energy.
  - 'Thermodynamics.txt'. Four column file with time (in time steps), kinetic energy, potential energy and instant temperature.
  - 'VelInicial.txt'. Contains the initial velocity components in one column.
  - 'VelFinal.txt'. The same as 'VelInicial.txt' but for final velocities.

The simulation flow is the following:
  - The system starts at a BCC crystalline structure which allows high density configurations.
  - It is disordered running a small simulation with a thermal bath (Andersen's thermostat with reduced temperature 100 and velocity Verlet integrator).
  - The velocities are initialize with a bimodal distribution.
  - The system is thermalized by performing a small simulation whitout measurements.
  - Production loop measuring after every integration step.

## Whit thermal Reservoir

In this section, I performed the simulation in the (T, V, N) ensemble so a thermal bath is implemented in the main section fo the code (Andersen's thermostat and velocity verlet integrator).

The production code is called 'Codigo.cpp'. For a correct run of the code, it requires three directories: 'Thermodynamics', 'Distribution', 'Difusion'. It takes a set of values of the density (in reduced units), the scale of energy, length and mass and the thermostat temperature and produces three files: 
  - 'Thermodynamics/Density = <The one seleced>.txt'. Five column file with time (in time steps), kinetic energy, potential energy,  temperature and pressure (in standard units). Used to measure the equilibrium properties of the fluid at a given temperature.
  - 'Difusion/Density = <The one seleced>.txt'. Two column file with time and average square displacement of the particles. Used to obtain the diffusion coefficent through a linear fit.
  - 'Distribution/Density = <The one seleced>.txt' One column file with the distance between particles. Used to obtain the radial distribution function.

The simulation flow is the following:
  - The system starts at a BCC crystalline structure which allows high density configurations.
  - It is disordered running a small simulation with a thermal bath (Andersen's thermostat with reduced temperature 100 and velocity Verlet integrator).
  - The velocities are initialize with a gaussian distribution with variance equal to thermal bath temperature.
  - The system is thermalized by performing a small simulation whithout measurements.
  - Production loop. Measuremets are done every 0.05 time steps (which correspond to 50 integration steps).

The analysis code is called 'Analisis_code.cpp'. It requires a directory called 'Analisis'. It uses the binning method to obtain the average values and its error of the potential energy, kinetic energy, temperature and pressure at equilibrium with a given density. It produces 4 files, one for each magnitude.
