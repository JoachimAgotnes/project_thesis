function [test_data, train_data] = getData(well_nr)
par = ParametersGasLift; 
%% Getting training set
% load('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\data_noise_train.mat');
disp('Choose a file for training data');
file = uigetfile('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\');
load(['C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\',file]);
y = z.OutputData;
u = z.InputData;
y_train = y(:,well_nr);
u_train = u(:,well_nr);
train_data = iddata(y_train,u_train,par.T);

%% Getting test data
% load('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\data_noise_test.mat');
disp('Choose a file for test data');
file = uigetfile('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\');
load(['C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\',file]);
y = z.OutputData;
u = z.InputData;
y_test = y(:,well_nr);
u_test = u(:,well_nr);
test_data = iddata(y_test,u_test,par.T);
end