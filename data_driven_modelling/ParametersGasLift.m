function  par = ParametersGasLift

%number of wells
par.n_w = 3;
%gas constant
par.R = 8.314; %[m3 Pa/(K mol)]
%molecular weigth
par.Mw = 20e-3; %[kg/mol] -- Attention: this unit is not usual

%% Properties
%density of oil - dim:  nwells x 1
par.rho_o = 8*1e2;%[8;8].*1e2; %[kg/m3] 
%riser oil density 
par.rho_ro = par.rho_o; %[kg/m3] 
%1cP oil viscosity
par.mu_oil = 1*0.001; %[Pa s or kg/(m s)] 

%% Project
%well parameters - dim:  nwells x 1
%length
par.L_w = [1500;1500;1500]; %[m]
%height
par.H_w  = [1000;1000;1000]; %[m]
%diameter
par.D_w = [0.121;0.121;0.121]; %[m]
%well transversal area
par.A_w = pi.*(par.D_w/2).^2;%[m2]

%well below injection - [m]
par.L_bh = [500;500;500];
par.H_bh = [500;500;500];
par.D_bh = [0.121;0.121;0.121];
par.A_bh = pi.*(par.D_bh/2).^2;%[m2]

%annulus - [m]
par.L_a = par.L_w;
par.H_a = par.H_w;
par.D_a = [0.189;0.189;0.189];
%volume of the annulus
par.V_a = par.L_a.*(pi.*(par.D_a/2).^2 - pi.*(par.D_w/2).^2); %[m3]

%riser - [m]
par.L_r = 500;
par.H_r = 500;
par.D_r = 0.121;
%riser areas
par.A_r = pi.*(par.D_r/2).^2;%[m2]

%injection valve characteristics  - dim:  nwells x 1
par.C_iv = [0.1e-3;0.1e-3;0.1e-3];%[m2]
%production valve characteristics  - dim:  nwells x 1
par.C_pc = [2e-3;2e-3;2e-3];%[m2]
%riser valve characteristics
par.C_pr = [10e-3];%[m2]
%parameter to account for differences in gas and liquid pressures
par.slip_real = 1;

%parameters
%reservoir pressure
par.p_res = [150;155;160]; % [bar] 
%Annulus temperature
par.T_a = [28+273;28+273;28+273]; %[K]
%well temperature
par.T_w = [32+273;32+273;32+273]; %[K]
%riser temperature
par.T_r = 30+273; %[K]
%separator pressure
par.p_s = 20; %[bar]

%sampling time /control interval /1 simulation iteration time
par.T = 3600*24; % [s]

%% For erosion model
% Sand
par.d_p = 2.5*10^(-4); %[m] particle diameter
par.rho_p = 2.5*10^3; %[kg/m3] particle density
par.mdot_p = 100*10^(-3); %[kg/s] sand rate

% Choke
par.K = 2*10^(-9); %[-] material erosion constant
par.rho_t = 7800; %[kg/m3] sensity CS
par.r = 0.2; %[m] radius of curvature
par.D = 0.05; %[m] Gap between body and cage
par.H = 0.15; %[m] Height of gallery

% Constants
par.C_unit = 1000; % Unit conversion factor: now in mm/s
par.C_1 = 1.25; %[-] Model/geometry factor
par.n = 2.6; %[-] Velocity coefficient
par.GF = 2; %[-] Geometry factor

% Precalculations of erosion in choke:
par.alpha = atan(1/sqrt(2*par.r));
par.F = 0.6*(sin(par.alpha) + 7.2*(sin(par.alpha) - sin(par.alpha)^2))^0.6 * (1-exp(-20*par.alpha));
par.A_g = 2*par.H*par.D; %[m2] Effective gallery area
par.A_t = par.D_w(1)^2*pi/(4*sin(par.alpha)); % Area exposed to erosion
par.ER_constant = par.K*par.F*par.C_1*par.GF*par.mdot_p *par.C_unit/(par.rho_t*par.A_t);


%System parameters for nominal model
par.GOR = [0.10;0.12;0.11];
par.PI = [5;5;5];


