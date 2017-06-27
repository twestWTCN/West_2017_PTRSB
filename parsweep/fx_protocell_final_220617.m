function [xdot] = fx_protocell_final_220617(x,Veq_crys_cyto,R_org_cat,K_aa,V_cyto,dt)
%%%%%%%%%%%%%%%%%%%% MAIN PROTOCELL FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function sets the fixed parameters and the system of ODEs that
% describe the dynamics of the protocell model.
%
% West et al. (2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Parameters %%%%%%%%%
% Constants
AN = 6.02e23;                   % mol^-1: Avagadro's constant
% Radii (dm)
r_mem = 1e-7;                   % dm: - Cell membrane thickness ~10nM
% Areas (dm^2)
% phi_fa = 6e-17;                % dm^2: - fatty acid head group size
phi_fa = 2e-17;                % dm^2: - fatty acid head group size 0.2 nm^2
% Volumes (dm^3)
V_crys_min = 1e-19;             % dm^3: - minimum crystal size

% Permeabilities (dm.s^-1)
P_crys_surf = 1e-13;            % dm.s^-1: - Crystal diffusion from bulk to crystal 
P_crys_cyto = 1e-13;            % dm.s^-1: - Crystal diffusion from cyto
P_crys_mem = P_crys_cyto;       % dm.s^-1: - Crystal diffusion from membrane
P_aa_cyto = 1e-10;              % dm.s^-1: - AA diffusion from cyto
P_fa_cyto = 1e-10;             % dm.s^-1: - fa diffusion from cyto

% Saturation Constants (mol.dm^-3)
K_CO2 = 3e-4;                   % mol.dm^-3: binding constant of CO2 for FeS
K_aa = K_aa;                    % mol.dm^-3: binding constant of AA for FeS
K_crys = 1e-8;                  % mol.dm^-3: saturation constant of FeS from bulk

% Rate Constants (s^-1)
k_grow = 1e-6;                  % s^-1: rate constant for evolution of crystal size distribution

% Concentrations (mol.dm^-3)
C_CO2_cyto = 1e-3;              % mol.dm^-3: Concentration of aqueous CO2 in cyto
C_aa_sink = 1e-6;               % mol.dm^-3: Concentration of aa in sink
C_fa_sink = 1e-6;              % mol.dm^-3: Concentration of fa in sink
C_crys_sink = 1e-9;             % mol.dm^-3: - concentration in sink

% Areal Turnover (mol.dm^-2.s^-1)
R_org_cat = R_org_cat;          % mol.dm^-2.s^-1: Areal molar turnover of organic production

% Organic Fraction (unitless)
lambda_cys = 0.1;               % fraction of cysteine produced
lambda_fa = 0.25;              % fraction fatty acid produced

%% %%%% Intialise Variables %%%%%%
V_crys_cyto = x(1);             % dm^3: Crystal population mean volume
C_crys_mem = x(2);              % mol.dm^-3: Concentration of crystals in the membrane
C_aa_cyto = x(3);               % mol.dm^-3: Concentration of amino acids in the cytoplasm
C_fa_cyto = x(4);              % mol.dm^-3: Concentration of fatty acids in the cytoplasm
SA_cyto = x(5);                 % dm^2: Surface area of the cytoplasm


% Stepwise recomputation of variables
% Amino Acid Modulator
eps_aa = ((C_aa_cyto^1)/((C_aa_cyto^1)+K_aa));      % Unitless - Saturating Michaelis Menten kinetics of amino acid binding FeS
% Cell Geometry
r_cyto = sqrt(SA_cyto/(4*pi));                      % dm: - Radius of spherical cell 
V_mem = (4/3)*pi*((r_cyto+r_mem)^3 - r_cyto^3);     % dm^3: volume of membrane shell

% Crystal Dynamics
N_crys_cyto = Veq_crys_cyto/V_crys_cyto;            % #: number of crystals in cytosol
C_crys_cyto = N_crys_cyto/(AN*V_cyto);              % mol.dm^-3: - concentration of crystal in cytosol
SA_crys = 6*(V_crys_cyto^(2/3));                    % dm^2: cuboidal surface area of crystals

%% %%%% Main Dynamics Begin Here %%%%%%%%%%
% Crystal Flows and Volume Changes
% Crystal transport
dC_crys_cyto_sink   =    P_crys_cyto*(SA_cyto/V_cyto)*(C_crys_sink-C_crys_cyto);               % mol.dm^-3.s^-1: - diffusion out of cell
dC_crys_mem_cyto    =    P_crys_cyto*(SA_cyto/V_mem)*(C_crys_cyto-C_crys_mem)*(eps_aa);         % mol.dm^-3.s^-1: - association of crystal with membrane
dC_crys_mem_sink    =    P_crys_cyto*(SA_cyto/V_mem)*(C_crys_sink-C_crys_mem)*(1-eps_aa);      % mol.dm^-3.s^-1: - dissociation of crystal membrane to sink
dC_crys_mem = dC_crys_mem_sink + dC_crys_mem_cyto; % total rate of change to membrane concentration

% Crystal deposition
dC_crys_cyto_growth = P_crys_surf*((N_crys_cyto*SA_crys)/V_cyto)*(C_crys_cyto-K_crys)*(1-eps_aa);   % mol.dm^-3.s^-1: - lumped growth/agglomeration  

% Find changes in crystal concentration with respect to the cytoplasm
c_gain = dC_crys_cyto_growth;   % gains are always from growth
c_loss = [];

if dC_crys_cyto_sink>0;     c_gain = c_gain; else    c_loss = -dC_crys_cyto_sink;        end % add negative flow (leak)
if dC_crys_mem_cyto<0 ;     c_gain = c_gain; else    c_loss = c_loss+dC_crys_mem_cyto;   end % add flow (membrane association)    
if isempty(c_loss);  c_loss = 1e-20; end % If all gains then make losses very small

mbal = c_gain/c_loss; % ratio of gains to losses

if (V_crys_cyto-V_crys_min)>1e-20; 
        dV_crys_cyto = k_grow*(mbal-1)*(V_crys_cyto-V_crys_min); % dm^3.s^-1: evolution of mean size of crystal population
else
        dV_crys_cyto = 0; % if difference very small then do not compute (ensures numerical stability)
end
    
% Catalysis and Cell Growth
% organic synthesis
R_org_cat    =   R_org_cat/AN;                              % convert from areal turnover (cm^-2.s^-1) to areal molar (mol.cm^-2.s^-1)
Ncrys_mem    =   C_crys_mem*V_mem*AN;                       % #: number of crystals in membrane
v_orgs       =   (SA_crys*Ncrys_mem*R_org_cat)/(V_cyto);    % mol.dm^-3.s^-1: maximum rate of organic synthesis
dC_orgs_cyto =   v_orgs*(C_CO2_cyto/(C_CO2_cyto+K_CO2));    % mol.dm^-3.s^-1: substrate dependent (CO2) rate of synthesis

% Individual organic partitioning + leak permeability
dC_cys_cyto     =   lambda_cys*dC_orgs_cyto + ((SA_cyto/V_cyto)*P_aa_cyto*(C_aa_sink-C_aa_cyto));       % mol.dm^-3.s^-1: rate of concentration change of aa
dC_fa_cyto     =   lambda_fa*dC_orgs_cyto + ((SA_cyto/V_cyto)*P_fa_cyto*(C_fa_sink-C_fa_cyto));    % mol.dm^-3.s^-1: rate of concentration change of fa

% Membrane Growth
dSA_cyto    =   (1/2)*dC_fa_cyto*V_cyto*AN*phi_fa;     % cm^2.s^-1: rate of change of membrane surface area

xdot(1) = dV_crys_cyto;
xdot(2) = dC_crys_mem;
xdot(3) = dC_cys_cyto;
xdot(4) = dC_fa_cyto;
xdot(5) = dSA_cyto;
xdot(6) = C_crys_cyto;
