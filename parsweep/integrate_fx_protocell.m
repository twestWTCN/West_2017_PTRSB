function xcept = integrate_fx_protocell(R_orgs_cat,K_AA,N,dt,div) 
%%%%%%%%%%%%%%%%%%%% SIMULATION SETUP AND INTEGRATOR %%%%%%%%%%%%%%%%%%%%%%
if nargin<5
    div = 0; % If unspecified - dont divide
end


%% Set Initial Conditions
% Cell geometry
V_cyto = 1e-12;                          % dm^3: V^cyto - 1e15 m^3 - 1000um^3
SA_cyto = (pi^(1/3))*((6*V_cyto)^(2/3)); % dm^2: SA^cyto- SA of V_cyto

% Initial States
x(1) = 1e-18;   % dm^3: V_crys_cyto - crystal of 1um length  - typical mackinawite microcrystal (Harmandas et al. 2006)
x(2) = 1e-9;    % mol.dm^-3: C_crys_mem -  (set low)
x(3) = 1e-6;    % mol.dm^-3: C_aa_cyto- mol.dm^-3 (set low)
x(4) = 1e-6;    % mol.dm^-3: C_lip_cyto - mol.dm^-3 (set low)
x(5) = SA_cyto; % dm^2: SAcell -  
x(6) = 1e-6;    % ambient crystal mol.dm^-3

% Equilibirum crystal pop volume
AN = 6.02e23;
Veq_crys_cyto = AN*(x(6)*V_cyto)*x(1) ;

xstore(:,1) = x;
xstore(5,1) = x(5)/100; % convert surface area (dm^2 to cm^2)

% Euler integration
for i = 2:N
    fx =  fx_protocell_final_220617(x,Veq_crys_cyto,R_orgs_cat,K_AA,V_cyto,dt);
    x(:,1:5) = x(:,1:5) +(dt*fx(:,1:5) ); % x + (dx/dt * dt)
    x(:,6) = fx(:,6);

    % If division
    if div == 1
        if x(5)> 1e-7
            % Reset Initial Conditions
            x(1) = x(1);
            x(2) = x(2)/2;
            x(3) = x(3)/2;
            x(4) = x(4)/10;
            x(5) = x(5)/2;
            x(6) = x(6);
            disp('Cell Divided!!!')
        end
    end
    xstore(:,i) = x;
    xstore(5,i) = x(5)/100; % convert surface area (dm^2 to cm^2)
% disp((i/N)*100)
end
 for i = 1:size(xstore,1)
     xcept(i,:) = downsample(xstore(i,:),10); % Downsample to reduce memory load
 end

