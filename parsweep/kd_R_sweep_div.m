%%%%%%%%%%%%%%%%%%%% Dividing Protocell Parameter Sweep %%%%%%%%%%%%%%%%%%%
% This script will run similations across a 2D matrix of parameters given % 
% by R_orgs_cat (areal turnover rate) and K_aa (amino acid binding        %
% constant. These simulations will allow the protocell to divide in the   %
% right conditions.Integrated series are saved for graphical output.      %
% For figure 2 of West et al. (2017)                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
t_simend = 200*24*60*60;
dt = 60;
N = t_simend/dt;
tvec =  linspace(0,t_simend/60/60/24,N);

Nmn = 16;
 
R_orgs_cat = 10.^linspace(11,14,Nmn);
K_aa = 10.^linspace(-5.5,-2.5,Nmn); % log10 space vector of binding constants

for n =1:Nmn
    x = zeros(6,1);
    parfor m = 1:Nmn
        xstore = integrate_fx_protocell(R_orgs_cat(n),K_aa(m),N,dt,1);
       MNxstore{m,n} = xstore;
       disp([m n])
    end
end
AN = 6.02e23;
R_orgs_cat = R_orgs_cat./AN; % convert to molar
tvec = downsample(tvec,10);

titname = {'Crystal size', 'Crystals in the membrane','Concentration of amino acids','Concentration of lipids','Protocell surface area','Crystals in cytosol'};

if ~exist([pwd '\parsweep\saves\']); mkdir([pwd '\parsweep\saves\']); end
save([ pwd '\parsweep\saves\divparlists'],'Nmn','R_orgs_cat','K_aa')
save([pwd '\parsweep\saves\divsimpar'],'t_simend','dt','N','tvec','titname')
savesweepparts('divided',Nmn,MNxstore)
