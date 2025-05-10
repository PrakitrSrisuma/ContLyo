% ==============================================================================
% This is a routine for desorption/diffusion analysis.
% Continuous lyophilization.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Pre-simulation
% Add paths
addpath('Input Data', 'Model Equations', 'Events','Exporting Graphics','Plotting', ...
    'Validation Data','Simulations','Calculations');

mode = 0;

if mode == 0
    l = 5e-7;
    Deff = 7e-16;
    teff = l^2/Deff;
    keff = 7.8e-5;
    teff2 = 1/keff;
    
    fprintf(['Time scale for diffusion is ' num2str(teff) '\n'])
    fprintf(['Time scale for desorption is ' num2str(teff2) '\n'])

elseif mode == 1
    
    l = 1e-2;

    eps = .815;
    tau = 1.2;
    M = 18;
    T = 256;
    D12 = 1.97e-5;
    re = 5e-6;

    D12eff = (eps/tau)*D12;

    % v1 = (8*8.314e3*T/(pi*M))^0.5;
    % Dk = (2/3)*(eps/tau)*v1*re;
    Dk = 97*re*sqrt(T/M);
    Dall = 1/(1/D12eff + 1/Dk);

    teff = l^2/Dall;
    keff = 7.8e-5;
    teff2 = 1/keff;
    
    fprintf(['Time scale for diffusion is ' num2str(teff) '\n'])
    fprintf(['Time scale for desorption is ' num2str(teff2) '\n'])


end