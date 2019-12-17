clear 
close all
clc

well_nr = 3;
disp(['Choose the result file for ARX for well number ',num2str(well_nr)]);
file = uigetfile('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\');
arx_ = load(['C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\',file]);
results_arx = arx_.results_arx;

disp(['Choose the result file for ARMAX for well number ',num2str(well_nr)]);
file = uigetfile('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\');
armax_ = load(['C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\',file]);
results_armax = armax_.results_armax;

disp(['Choose the result file for OE for well number ',num2str(well_nr)]);
file = uigetfile('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\');
oe_ = load(['C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\',file]);
results_oe = oe_.results_oe;

par = ParametersGasLift;
[train_data,test_data] = getData(well_nr);

%% Finding the best model from aic. 
% ARX
[min_arx, id_arx] = min(results_arx(:,5));
[min_arx_fpe, id_arx_fpe] = min(results_arx(:,6));
if id_arx ~= id_arx_fpe
    disp('AIC and FPE for ARX does not give the same best model');
end  
best_arx = results_arx(id_arx,:);
sys_arx = arx(train_data,best_arx(2:4));
subplot(3,1,1);
compare(test_data,sys_arx); 
title('');
legend('Test data','ARX model');

% ARMAX
[min_arm, id_arm] = min(results_armax(:,6));
[min_arm_fpe, id_arm_fpe] = min(results_armax(:,7));
if id_arm ~= id_arm_fpe
    disp('AIC and FPE for ARMAX does not give the same best model');
end  
best_arm = results_armax(id_arm_fpe,:);
sys_arm = armax(train_data, best_arm(2:5));
subplot(3,1,2);
compare(test_data,sys_arm);
title('');
legend('Test data','ARMAX model');


% OE
[min_oe, id_oe] = min(results_oe(:,5));
[min_oe_fpe, id_oe_fpe] = min(results_oe(:,6));
if id_oe ~= id_oe_fpe
    disp('AIC and FPE for OE does not give the same best model');
end  
best_oe = results_oe(id_oe,:);
sys_oe = oe(train_data, best_oe(2:4));
subplot(3,1,3);
compare(test_data,sys_oe);
title('');
legend('Test data','OE model');

figure('Name','ARX');
resid(train_data,sys_arx);
ylabel(['Predicted erosion in well ',num2str(well_nr)]);

fig2 = figure('Name','ARMAX');
resid(train_data,sys_arm);
ylabel(['Predicted erosion in well ',num2str(well_nr)]);

figure('Name','OE');
resid(train_data,sys_oe);
ylabel(['Predicted erosion in well ',num2str(well_nr)]);


