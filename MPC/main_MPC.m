clear all
close all
clc 
addpath("C:\Users\Joachim\Documents\casadi-windows-matlabR2016a-v3.4.5")
import casadi.*

dt = 3600; 
time_total = 3600*24*365; % seconds
tot_it = 500;
[x0,z0,u0] = InitialConditionGasLift_5;
x_mes_next = x0;
z_mes_next = z0;
x_next = x0;
z_next = z0;
u_k = u0;
u_opt = [];
par = ParametersGasLift;
x_vec = [];
obj_vec = [];
z_vec = [];
tic;
nm = 70;
np = 100;
f = waitbar(0,'MPC'); % For visual display of progress
simLength = 500;

s_vec = [];
w_vec = [];
mu = 0;
sigma = 1;

flag = [];



for t = 0:simLength
   waitbar(t/tot_it,f, [num2str(round(t/tot_it*100,2)),'%... Iteration ', num2str(t),'...    Elapsed time: '...
   ,datestr(seconds(toc),'HH:MM:SS')]);
    % Calculating input with NMPC 
    [u_k,s,w_,exitflag]=  NMPC(x_next, u_k,z_next,nm,np);
    flag(t+1) = exitflag;
    % Adding optimal input to vector
    s_vec = [s_vec,s];
    w_vec = [w_vec,w_];
    u_opt = [u_opt, u_k];

    % Finding the state after applying input u
    [x_next,z_next] = WellPlantModel(x_next,z_next,u_k,par);
    
    % Adding noice, if active, use these in NMPC
%     x_mes_next = x_next  + normrnd(mu,sigma,[3,1]).*x0*0.015;
%     z_mes_next = z_next + normrnd(mu,sigma,[62,1]).*z0*0.015;
    x_vec = [x_vec,x_next];
    obj_vec(t+1) = sum(z_next(28:30));
    z_vec = [z_vec,z_next];
end
close(f)

%% Plotting
t = [0:1:simLength];

subplot(3,1,1);
stairs(t,u_opt(1,:));
hold on
stairs(t,u_opt(2,:));
stairs(t,u_opt(3,:));
legend('Well 1','Well 2','Well 3');
ylim([0,2.2]);
xlabel('Time [day]');
ylabel('Gas lift rate [kg/s]');

subplot(3,1,2);
plot(t,x_vec(1,:));
hold on
plot(t,x_vec(2,:));
plot(t,x_vec(3,:));
legend('Well 1','Well 2','Well 3');
xlabel('Time [day]');
ylabel('Erosion [mm]');

subplot(3,1,3);
plot(t,obj_vec);
xlabel('Time [day]');
ylabel('Total production of oil [kg/s]');

