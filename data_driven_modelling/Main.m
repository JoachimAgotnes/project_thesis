clear  
close all
clc
%% Simulation
%paramters
par = ParametersGasLift;
simLength = 500;
%initial condition
[dx0,z0,u0] = InitialConditionGasLift_5;

%System parameters for nominal model
GOR = par.GOR;
PI = par.PI;

%states to measurement mapping function
H = zeros(39,62);
H(1,1) = 1; %pa1 - annulus pressure, well 1
H(2,2) = 1; %pa2 - annulus pressure, well 2
H(3,3) = 1; %pa3 - annulues pressure, well 3
H(4,4) = 1; %pw1 - well head pressure, well 1
H(5,5) = 1; %pw2 - well head pressure, well 2
H(6,6) = 1; %pw2 - well head pressure, well 3
H(7,7) = 1; %pw1 - injection point pressure, well 1
H(8,8) = 1; %pw2 - injection point pressure, well 2
H(9,9) = 1; %pw3 - injection point pressure, well 3
H(10,10) = 1; %pbh1 - bottom hole pressure, well 1
H(11,11) = 1; %pbh1 - bottom hole pressure, well 2
H(12,12) = 1; %pbh1 - bottom hole pressure, well 3
H(13,16) = 1; % mixture denisty in tubing 1
H(14,17) = 1; % mixture denisty in tubing 2
H(15,18) = 1; % mixture denisty in tubing 3
H(16,37) = 1; %p_rh - riser head pressure
H(17,39) = 1; %p_m - manifold pressure
H(18,41) = 1; %w_to - total oil production
H(19,42) = 1; %w_tg - total gas production
H(20,43) = 1; %gas holdup @ annulus 1
H(21,44) = 1; %gas holdup @ annulus 2
H(22,45) = 1; %gas holdup @ annulus 3
H(23,46) = 1; %gas holdup @ well 1
H(24,47) = 1; %gas holdup @ well 2
H(25,48) = 1; %gas holdup @ well 3
H(26,49) = 1; %oil holdup @ well 1
H(27,50) = 1; %oil holdup @ well 2
H(28,51) = 1; %oil holdup @ well 3
H(29,52) = 1; %gas holdup @ riser
H(30,53) = 1; %oil holdup @ riser
H(31,54) = 1; %V_p @ well 1
H(32,55) = 1; %V_p @ well 2
H(33,56) = 1; %V_p @ well 3
H(34,57) = 1; %mixed viscosity @ well 1
H(35,58) = 1; %mixed viscosity @ well 2
H(36,59) = 1; %mixed viscosity @ well 3
H(37,60) = 1;
H(38,61) = 1;
H(39,62) = 1;
ny = 39; %number of measurements

%% Initialization
counter = 1; %counter for the number of executions

%Simulation
dxk = dx0;
zk = z0;
uk = u0;

%For plotting
inputPlantArray = u0;
yPlantArray = H*z0;
ER_array = [0;0;0];
rho_p = par.rho_p; %[kg/m3] particle density
r = par.r;
alpha = par.alpha;
d_p = par.d_p; %[m] particle diameter
gamma_ = d_p./par.D_w;

%% Exciting the system

% For trainng
uk = repmat(u0,1,simLength);
uk(1:end) = uk(1:end)*4;
% uk(:,201:350) = uk(:,201:350) + 0.5*uk(1);
% uk(:,351:end) = uk(:,351:end) - 0.5*uk(1);

% For test set
% uk = repmat(u0,1,simLength);
% uk(:,1:150) = uk(:,1:150) - 0.3*uk(1);
% uk(:,201:300) = uk(:,201:300) + 1*uk(1);

%% Simulation
gamma_counter = 0;
A = [];
gamma_c = [];
sensitivities = [];

mu = 0;
sigma = 0.01;
while counter <= simLength
  fprintf('>>> Iteration: %d \n',counter)  
  
%1. simulate plant (SS) and store plant data
    [dxk,zk] = WellPlantModel(dxk,zk,uk(1:3,counter),par);
    
    % Adding measurement noise 
    dxk1 = dxk + normrnd(mu,sigma,[3,1]);
    
    for i = 1:3
        if dxk1(i) < 0
            dxk1(i) = 0;
        end
    end
     ER_array = [ER_array, dxk1];  
        %%%%%%%%%%%%%%%%% plant information
        inputPlantArray = [inputPlantArray,uk(1:3,counter)];
        yPlantArray = [yPlantArray,H*zk];
        %%%%%%%%%%%%%%%%%%
   
        rho_m = yPlantArray(13:15,counter+1);
        V_p = yPlantArray(31:33,counter+1);
        mu_m = yPlantArray(34:36,counter+1);
        A(1:3,counter) = rho_m.^2*tan(alpha).*V_p.*par.D./(rho_p.*mu_m);
        gamma_c(1:3,counter) = rho_m./(rho_p.*(1.88.*log(A(1:3,counter))-6.04));
        
        % Checking if any gammas are over 0
        if gamma_c(1,counter) > 0 || gamma_c(2,counter) > 0 || gamma_c(3,counter) > 0
            disp(['Gamma over 0 in iteration: ',counter]);
            gamma_counter = gamma_counter + 1;
        end

        counter = counter + 1;
        
end

disp(['Number of positive gammas: ',num2str(gamma_counter)]);

%% Saving the results
time = 0:1:simLength; %[s]
plot(time,ER_array(1,:),'LineWidth',1);
hold on
plot(time,ER_array(2,:),'LineWidth',1);
plot(time,ER_array(3,:),'LineWidth',1);
xlabel('Time [days]');
ylabel('Erosion [mm]');
legend('Well 1','Well 2','Well 3');
z = iddata(ER_array(:,2:end)',uk',par.T);
% save('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\data_train.mat','z');







