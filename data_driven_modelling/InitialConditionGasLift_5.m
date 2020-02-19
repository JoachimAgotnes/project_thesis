function [dx0,z0,u0] = InitialConditionGasLift_5

%% Differential states

%well erosion rate 
ER0 = [1,1,1]'/(365*24*3600); %[mm/s] 9-11

%% Algebraic states
%pressure - annulus
p_ai0 = [61.9230, 62.1454, 62]';%[bar] 1-3  (bar to Pa  = x10^5)
%pressure - well head
p_wh0 = [42.5851, 45.3082, 44]';%[bar] 4-6
%pressure - injection point
p_wi0 = [56.8713, 57.1119, 57]';%[bar] 7-9
%pressure - below injection point (bottom hole)
p_bh0 = [96.1433, 98.3867, 96.25]';%[bar] 10-12
%density - annulus
rho_ai0 = [0.4949, 0.4967, 0.4955]';%[100 kg/m3] 13-15 
%mixture density in tubing
rho_m0 = [2.3400, 2.2359, 2.3]';%[100 kg/m3] 16-18
%well injection flow rate
w_iv0 = [0.5000, 0.5000,0.5000]';%[kg/s] 19-21 
%wellhead total production rate
w_pc0 = [30.1212, 33.3235,32]';%[kg/s] 22-24 
%wellhead gas production rate
w_pg0 = [3.1928, 4.0168, 4]';%[kg/s] 25-27 
%wellhead oil production rate
w_po0 = [26.9283, 29.3067, 28]';%[kg/s] 28-30 
%oil rate from reservoir
w_ro0 = [26.9283, 29.3067, 28]';%[kg/s] 31-33 
%gas rate from reservoir
w_rg0 = [26.9283, 35.1680, 28]';%[0.1 kg/s] 34-36 
%riser head pressure
p_rh0 = 22.9558;%[bar] 37 
%mixture density in riser
rho_r0 = 1.3618;%[100 kg/m3] 38 
%manifold pressure
p_m0 = 32.8920;%[bar] 39
%riser head total production rate
w_pr0 = 63.4446;%[kg/s] 40 
%riser head total oil production rate
w_to0 = 56.2350;%[kg/s] 41 
%riser head total gas production rate
w_tg0 = 7.2096;%[kg/s] 42 

%%setting diff states as algebraic
%gas holdup @ annulus
m_ga0 = [1.0568, 1.0606, 1.0644]';%[ton] 43-45 mga(2) +(mga(2)-mga(1))
%gas holdup @ well
m_gt0 = [0.7470, 0.7956, 0.8442]';%[ton] 46-48 mgt(2) + (mgt(2)-mgt(1))
%oil holdup @ well
m_ot0 = [6.3000, 5.8047, 5.3094]';%[ton] 49-51
%gas holdup @ riser
m_gr0 = 0.1265;%[ton] 52  
%oil holdup @ riser
m_or0 = 0.9863;%[ton] 53
%particle impact velocity
V_p0 = [2,2,2]'; % 54-56
%mixed dynamic viscosity
mu_f0 = [0.001,0.001,0.001]'; % 57-59
%g1 
g10 = [0.2,0.2,0.2]'; % 60-62

%% Inputs
%gas lift rate
w_gl0 = [0.5,0.5,0.5]'; %[kg/s]

dx0 = vertcat(ER0);
z0 = vertcat(p_ai0,p_wh0,p_wi0,p_bh0,rho_ai0,rho_m0,w_iv0,w_pc0,w_pg0,w_po0,...
    w_ro0,w_rg0,p_rh0,rho_r0,p_m0,w_pr0,w_to0,w_tg0,m_ga0,m_gt0,m_ot0,m_gr0,m_or0,V_p0,mu_f0,g10);
u0 = w_gl0;