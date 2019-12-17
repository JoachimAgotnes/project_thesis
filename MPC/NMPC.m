function [u_,s,w_,exitflag] = NMPC(x_current,ui,z_current,nm,np)

addpath("C:\Users\Joachim\Documents\casadi-windows-matlabR2016a-v3.4.5")
import casadi.*
nu = 3; 
nx = 3; 
nz = 62; 

%% Parameters
par = ParametersGasLift;
%number of wells
n_w = 3; %[]
%gas constant
R = par.R; %[m3 Pa /K /mol]
%molecular weigth
Mw = par.Mw; %[kg/mol]

%properties
%density of oil - dim:  nwells x 1
rho_o = par.rho_o; %[kg/m3]
%riser oil density 
rho_ro = par.rho_ro;%[kg/m3]
%1cP oil viscosity
mu_oil = par.mu_oil;% [Pa s] 

%project
%well parameters - dim:  nwells x 1
L_w = par.L_w; %[m]
H_w = par.H_w; %[m]
D_w = par.D_w; %[m]
A_w = par.A_w;%[m2]

%well below injection - [m]
L_bh = par.L_bh;
H_bh = par.H_bh;
D_bh = par.D_bh;
A_bh = par.A_bh;%[m2]

%annulus - [m]
H_a = par.H_a;
V_a = par.V_a; %[m3]

%riser - [m]
L_r = par.L_r;
H_r = par.H_r;
D_r = par.D_r;
A_r = par.A_r;%[m2]

%injection valve characteristics  - dim:  nwells x 1
C_iv = par.C_iv;%[m2]
%production valve characteristics  - dim:  nwells x 1
C_pc = par.C_pc;%[m2]
%riser valve characteristics
C_pr = par.C_pr;%[m2]
% account for differences in the vapor and oil velocity
slip = par.slip_real;


%% For erosion model
% Sand
d_p = par.d_p; %[m] particle diameter
rho_p = par.rho_p; %[kg/m3] particle density
mdot_p = par.mdot_p; %[kg/s] sand rate

% Choke
K = par.K; %[-] material erosion constant
rho_t = par.rho_t; %[kg/m3] sensity CS
r = par.r; %[m] radius of curvature
D = par.D; %[m] Gap between body and cage
H = par.H; %[m] Height of gallery

% Constants
C_unit = par.C_unit; % Unit conversion factor: now in mm/s
C_1 = par.C_1; %[-] Model/geometry factor
n = par.n; %[-] Velocity coefficient
GF = par.GF; %[-] Geometry factor

% Precalculations of erosion in choke:
alpha = par.alpha;
F = par.F;
A_g = par.A_g; %[m2] Effective gallery area
ER_constant = par.ER_constant;
gma = d_p./D;

%% Algebraic states
%pressure - annulus
p_ai = MX.sym('p_ai',n_w);      % 1-3 [bar] (bar to Pa  = x10^5)
%pressure - well head
p_wh = MX.sym('p_wh',n_w);      % 4-6 [bar]
%pressure - injection point
p_wi = MX.sym('p_wi',n_w);      % 7-9 [bar]
%pressure - below injection point (bottom hole)
p_bh = MX.sym('p_bh',n_w);      % 10-12 [bar]
%density - annulus
rho_ai = MX.sym('rho_ai',n_w);  % 13-15 [100 kg/m3]
%mixture density in tubing
rho_m = MX.sym('rho_m',n_w);    % 16-18 [100 kg/m3]
%well injection flow rate
w_iv = MX.sym('w_iv',n_w);      % 19-21 [kg/s]
%wellhead total production rate
w_pc = MX.sym('w_pc',n_w);      % 22-24 [kg/s]
%wellhead gas production rate
w_pg = MX.sym('w_pg',n_w);      % 25-27 [kg/s]
%wellhead oil production rate
w_po = MX.sym('w_po',n_w);      % 28-30 [kg/s]
%oil rate from reservoir
w_ro = MX.sym('w_ro',n_w);      % 31-33 [kg/s]
%gas rate from reservoir
w_rg = MX.sym('w_rg',n_w);      % 34-36 [0.1 kg/s]
%riser head pressure
p_rh = MX.sym('p_rh',1);        % 37 [bar]
%mixture density in riser
rho_r =  MX.sym('rho_r',1);     % 38 [100 kg/s]
%manifold pressure
p_m = MX.sym('p_m',1);          % 39 [bar]
%riser head total production rate
w_pr = MX.sym('w_pr',1);        % 40 [kg/s]
%riser head total oil production rate
w_to = MX.sym('w_to',1);        % 41 [kg/s]
%riser head total gas production rate
w_tg = MX.sym('w_tg',1);        % 42 [kg/s]
%gas holdup @ annulus 
m_ga = MX.sym('m_ga',n_w); % 43-45 [ton]
%gas holdup @ well
m_gt = MX.sym('m_gt',n_w); % 46-48 [ton]
%oil holdup @ well 
m_ot = MX.sym('m_ot',n_w); % 49-51 [ton]
%gas holdup @ riser
m_gr = MX.sym('m_gr',1);   % 52 [ton]
%oil holdup @ riser
m_or = MX.sym('m_or',1);   % 53 [ton]
%particle impact velocity
V_p = MX.sym('V_p',n_w); % 54-56
%dynamic viscosity of mixture
mu_f = MX.sym('mu_f',n_w); % 57-59
g1 = MX.sym('g1',n_w); % 60-62


%control input
%gas lift rate
w_gl = MX.sym('w_gl',n_w); %[kg/s]

%parameters
p_res = MX.sym('p_res',n_w);
%productivity index
PI = MX.sym('PI',n_w); %[kg s^-1 bar-1]
%GasOil ratio
GOR = MX.sym('GOR',n_w); %[kg/kg]
%Annulus temperature
T_a = MX.sym('T_a',n_w); %[oC]
%well temperature
T_w = MX.sym('T_w',n_w); %[oC]
%riser temperature
T_r = MX.sym('T_r',1); %[oC]
%separator pressure
p_s = MX.sym('p_s',1); %[bar]
%time transformation: CASADI integrates always from 0 to 1 and the USER does the time
%scaling with T.
T = MX.sym('T',1); %[s]
%erosion rate
ER = MX.sym('ER',n_w); %

%% Modeling
%gas fraction (mass) of the well holdup - avoiding zero division
xGwH = (m_gt.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
%gas fraction (mass) of the riser holdup
xGrH = (m_gr.*1e3./(m_gr.*1e3+m_or.*1e3));

xGw = slip.*xGwH./(1 + (slip-1).*xGwH);
xOw = 1 - xGw;
xGr = slip.*xGrH./(1 + (slip-1).*xGrH);
xOr = 1 - xGr;

% ===================================
%     Well model with/withou pressure loss
% ===================================
% algebraic equations (all symbolic)
%annulus pressure - %g = 9.81
f1 = -p_ai.*1e5 + ((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*T_a).*((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a;
%well head pressure
f2 = -p_wh.*1e5 + ((R.*T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2;
%well injection point pressure
f3 = -p_wi.*1e5 + (p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w) + (128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3)));
%bottom hole pressure
f4 = -p_bh.*1e5 + (p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_ro./(3.14.*D_bh.^4.*rho_o));
%gas density in annulus
f5 = -rho_ai.*1e2 +(Mw./(R.*T_a).*p_ai.*1e5);
%fluid mixture density in well
f6 = -rho_m.*1e2 + ((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3);
%well injection flow rate
f7 = -w_iv + C_iv.*sqrt(rho_ai.*1e2.*max(0,(p_ai.*1e5 - p_wi.*1e5)));
%wellhead prodution rate
f8 = -w_pc + 1.*C_pc.*sqrt(rho_m.*1e2.*max(0,(p_wh.*1e5 - p_m.*1e5)));
%wellhead gas production rate
f9 = -w_pg + xGw.*w_pc;
%wellhead oil prodution rate
f10 = -w_po + xOw.*w_pc;
%oil from reservoir flowrate
f11 = -w_ro + PI.*1e-6.*(p_res.*1e5 - p_bh.*1e5);
%gas from reservoir production rate
f12 = -w_rg.*1e-1 + GOR.*w_ro;  
%riser head pressure
f13 = -p_rh.*1e5 + ((R.*T_r./Mw).*(m_gr.*1e3./(L_r.*A_r))) - ((m_gr.*1e3+m_or.*1e3 )./(L_r.*A_r)).*9.81.*H_r/2;
%riser density
f14 = -rho_r.*1e2 + ((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3);
%manifold pressure
f15 = -p_m.*1e5 + (p_rh.*1e5 + 9.81./(A_r.*L_r).*(m_or.*1e3+m_gr.*1e3).*H_r) + (128.*mu_oil.*L_r.*w_pr./(3.14.*D_r.^4.*((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3)));
%total production rate of well
f16 = -w_pr + 1.*C_pr.*sqrt(rho_r.*1e2.*(p_rh.*1e5 - p_s.*1e5));
%oil total production rate
f17 = -w_to + xOr.*w_pr; 
%gas total production rate
f18 = -w_tg + xGr.*w_pr; 
% setting differential equations as algebraic equations since the
% dynamics of ER is on a much larger time scale
f19 = (w_gl - w_iv).*1e-3;
f20 = (w_iv + w_rg.*1e-1 - w_pg).*1e-3;
f21 = (w_ro - w_po).*1e-3;
f22 = (sum(w_pg) - w_tg).*1e-3 ;
f23 = (sum(w_po) - w_to).*1e-3 ;
f24 = - V_p + 3/(4*A_g).*(w_po./rho_o + R*T_w.*w_pg./(p_wh.*10^5*Mw));
f25 = - mu_f + mu_oil.*(w_po./rho_o)./(w_po./rho_o + R.*T_w.*w_pg./(p_wh.*10^5*Mw));
f26 = -g1 + gma/0.1;
% differential equations - (all symbolic) - [ton]
% Erosion rate 
df1 =  ER_constant.*g1.*(V_p).^n; 

% Form the DAE system
diff = vertcat(df1);
alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26);

% give parameter values
alg = substitute(alg,p_res,par.p_res);
alg = substitute(alg,p_s,par.p_s);
alg = substitute(alg,T_a,par.T_a);
alg = substitute(alg,T_w,par.T_w);
alg = substitute(alg,T_r,par.T_r);

diff = substitute(diff,p_res,par.p_res);
diff = substitute(diff,T_w,par.T_w);

% concatenate the differential and algebraic states
x_var = vertcat(ER);
z_var = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg,p_rh,rho_r,p_m,...
    w_pr,w_to,w_tg,m_ga,m_gt,m_ot,m_gr,m_or,V_p, mu_f,g1);
p_var = vertcat(w_gl,GOR,PI,T);

%% NMPC 
umax = 2;
umin = 0.4;
dumax = 0.01;
x_threshold = 2; 


%% Defining empty nlp-problem
J = 0;
w = {};
w0 = [];
lbw =[];
ubw = [];
lbg = [];
ubg = [];
g = {};

x_prev = MX.sym('X0',nx);
Z0 = MX.sym('Z0',nz);
w = {w{:},x_prev}; % 1-3
lbw = [lbw,zeros(nx,1)]; 
ubw = [ubw,inf*ones(nx,1)];
w0 = [w0;x_current];

g = {g{:},x_prev};
lbg = [lbg,x_current];
ubg = [ubg,x_current];


%% Lifting initial conditions
uk = MX.sym('uk_init',nu);
w = {w{:}, uk}; % 4-6
w0 = [w0;  0.3*ones(nu,1)];
lbw = [lbw;ui];
ubw = [ubw;ui];


% uprev= MX.sym('uprev_init',nu);
U0 = MX.sym('U0',nu);
p = MX.sym('p',7);

x = MX.sym('x',nx);
uprev1 = uk;
uk1 = uk;
rho = 999999;
s = MX.sym('s',nx);
R = [1,0,0
    0,1,0
    0,0,1];
L =  -sum(w_po)+ rho*sum(s)  + 1/2 * ((uk1-U0)'*R*(uk1-U0));
f = Function('f',{x_var,z_var,p_var,uk1,U0,s},{diff,alg,L});

%% Using 3 collocation points:

h = par.T;
% Radau
t = [collocation_points(3, 'radau')];

% Finding M 
M = [t',0.5*t'.^2,1/3*t'.^3]*inv([[1;1;1],t',t'.^2]);

%% Looping through until timeend
for k = 1:np
    uprev = uk; 
    uk = MX.sym(['uk_' num2str(k)],nu);  
    w = {w{:}, uk}; % 7-9
    w0 = [w0;  ui];
    lbw = [lbw;umin*ones(nu,1)];
    ubw = [ubw;umax*ones(nu,1)];
    
    s = MX.sym(['s_',num2str(k)],nx);
    w = {w{:},s}; % 10-12
    lbw = [lbw;0*ones(nu,1)];
    ubw = [ubw;1*ones(nu,1)];
    w0 = [w0;0*ones(nu,1)];
    
    % Adding constraint for delta_u
    duk = uk - uprev;
    g = {g{:},duk};

    if k > nm
        lbg = [lbg;zeros(nu,1)];
        ubg = [ubg; zeros(nu,1)];
    else
        lbg = [lbg;-dumax*ones(nu,1)];
        ubg = [ubg;dumax*ones(nu,1)];
    end
    
    % Collocation points
    fk = [];
    Xk1 = [];
    zk = [];
    L1 = [];
    for d = 1:3
        Xk = MX.sym(['Xk_' num2str(k),'_',num2str(d)],nx);
        Zk = MX.sym(['Zk_' num2str(k),'_',num2str(d)],nz);
        Xk1 = [Xk1,Xk];
        w = {w{:}, Xk, Zk}; % 13-15
        w0 = [w0;x_current;z_current];
        lbw = [lbw;zeros(nx,1);0*ones(nz,1)];
        ubw = [ubw;inf*ones(nx,1);inf*ones(nz,1)];
     
        % Calculating xdot and objective function
        [fk1,zk1,L] = f(Xk,Zk,vertcat(uk,p),uk,uprev,s);
        L1 = [L1;L];
        if k > nm
            L1(d) = L - 1/2 * ((uk-uprev)'*R*(uk-uprev));
        end
        fk = [fk, fk1];
        zk = [zk,zk1];
    end
    ML = M*L1;
    J = J + ML(3);
%         ML = M(3,1)*sum(zk(28:30,1)) + M(3,2)*sum(zk(28:30,2)) + M(3,3)*sum(zk(28:30,3)) +rho*sum(s);
%         J = J + ML;
    x_next1 = [];
    for d = 1:3 
        % Calculating M*xdot for each collocation point
        Mfk = M(d,1)*fk(:,1) + M(d,2)*fk(:,2) + M(d,3)*fk(:,3);
        % Calculating x
        x_next = x_prev+h*Mfk;
        x_next1 = [x_next1,x_next];
        % Adding xk and Xk1 as constrains as they must be equal.
        g = {g{:},x_next-Xk1(:,d),zk(:,d)};
        lbg = [lbg;zeros(nx,1);zeros(nz,1)];
        ubg = [ubg;zeros(nx,1);zeros(nz,1)];  
    end

    % New NLP variable for state at end
    x_prev = MX.sym(['x_init_' num2str(k)],nx); 
    w = {w{:}, x_prev}; 
    w0 = [w0;x_current];
    lbw = [lbw;zeros(nx,1)];
    ubw = [ubw;inf*ones(nx,1)];
    
    % Gap
    g = {g{:},x_next-x_prev};
    lbg = [lbg;zeros(nx,1)];
    ubg = [ubg;zeros(nx,1)];
    
    % Constraint on erosion
    g = {g{:},x_prev-s};
    lbg = [lbg;zeros(nx,1)];
    ubg = [ubg;x_threshold*ones(nx,1)];
    
end

%% Solving optimization problem

% Formalizing problem 
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',vertcat(U0,p));

% Assigning solver (IPOPT)
solver = nlpsol('solver','ipopt',nlp);

% Solving the problem
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',[ui;par.GOR;par.PI;par.T]);

%% Extracting solution
w_opt = full(sol.x);
obj_opt = -full(sol.f); 
c_opt = full(sol.g);
lag_opt = full(sol.lam_g);
u_ = w_opt(7:9);
s = w_opt(10:12);
w_ = w_opt(43:45);
exitflag = solver.stats.success;
end