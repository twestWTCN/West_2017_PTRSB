%%% PROTOCELL MODEL MASTER SCRIPT %%%
addpath('parsweep')
% Following scripts are wrappers that call the main system of ODE's and Euler integration in
% 'fx_protocell_final.m' and 'integrate_fx_protocell.m' respectively

%% Simulations
kd_R_sweep                       % Run long time simulations to get equilibrium for range of parameters
kd_R_sweep_div               % Run simulations with divisions

%% Plotting
% Figure 2
waterfallplot_query % 2A
parspace_map % 2B & 2C

% Figure 3
waterfallplot_2D_div % 3
