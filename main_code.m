%% main_code.m
%%% APRIL 10, 2021

close all

%% Create the object 's' from the class 'simulator'
s = simulator;


%% Execute the simulation
s = s.simulate_model;


%% Plot the simulation
s.plot_simulation;

