%%%%%%%%%%%%%%%%%%%% Protocell Parameter Sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This script will run similations across a 2D matrix of parameters given % 
% by R_orgs_cat (areal turnover rate) and K_aa (amino acid binding        %
% constant. Integrated series are saved for graphical output              %
% For figures 1 of West et al. (2017)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all

% Set simulation parameters
t_simend = 120*(24*60*60);              % Simulation run time (days)(converted to seconds)
dt = 60;                                % Step size (s)
N = t_simend/dt;                        % Number of integration steps
tvec =  linspace(0,t_simend/60/60/24,N);% time vector (days)

% Define 2D space
Nmn = 16; % Resolution of parameter sweep 
R_orgs_cat = 10.^linspace(11,14,Nmn); % log10 space vector of turnover rates
K_aa = 10.^linspace(-5.5,-2.5,Nmn);   % log10 space vector of binding constants

%% Loop through simulations
for n = 1:Nmn
    parfor m = 1:Nmn
        xstore = integrate_fx_protocell(R_orgs_cat(n),K_aa(m),N,dt);
        MNxstore{m,n} = xstore;
        disp([m n])
    end
end

tvec = downsample(tvec,10); % Downsample to same as sample rate as integrated 

% Convert to molar turnover rate
AN = 6.02e23;
R_orgs_cat = R_orgs_cat./AN; % convert to molar


titname = {'Crystal size', 'Crystals in the protomembrane','Concentration of amino acids','Concentration of lipids','Protocell surface area','Crystals in cytosol'};
if ~exist([pwd '\parsweep\saves\']); mkdir([pwd '\parsweep\saves\']); end
   
save([pwd '\parsweep\saves\parlists'],'Nmn','R_orgs_cat','K_aa')
save([pwd '\parsweep\saves\simpar'],'dt','N','tvec','titname')

savesweepparts('query',Nmn,MNxstore)

