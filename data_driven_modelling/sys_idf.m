clear all
clc
close all
par = ParametersGasLift;
well_nr = 3;
%% Getting training and test data
[train_data, test_data] = getData(well_nr);

%% ARMAX and ARX models:
f = waitbar(0,'ARX and ARMAX'); % For visual display of progress
tot_it = 20*20*2*20; % total iterations 
curr_it = 0;
results_arx = []; % Matrix for storing results and parameters for ARX-models
results_armax = []; % Matrix for storing results and parameters for ARMAX-models
tic
for na = 1:1:20
    for nb = 1:1:20
        for nk = 0:1
           sys_arx = arx(train_data,[na,nb,nk]);
           [y_pred, fit, x0] = compare(test_data,sys_arx);
           results_arx = [results_arx; [fit,na,nb,nk, aic(sys_arx),fpe(sys_arx)]];
           
           for nc = 1:1:20
               sys_arm = armax(train_data,[na,nb,nc,nk]);
               [y_pred, fit, x0] = compare(test_data,sys_arm);
               results_armax = [results_armax; [fit,na,nb,nc,nk, aic(sys_arm),fpe(sys_arm)]];
               curr_it = curr_it+1;
               waitbar(curr_it/tot_it,f, ['ARX/ARMAX: ',num2str(round(curr_it/tot_it*100,2)),'%...    Elapsed time: '...
                   ,datestr(seconds(toc),'HH:MM:SS')]);
           end
        end
    end
end
close(f);
save('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\results_arx3.mat','results_arx');
save('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\results_armax3.mat','results_armax');

%% OE model
tic
results_oe = [];
f = waitbar(0,'OE');
tot_it = 20*20*2;
curr_it = 0;
for nb = 1:1:20
    for nf = 1:1:20
        for nk = 0:1
            sys_oe = oe(train_data,[nb,nf,nk]);
            [y_pred, fit, x0] = compare(test_data,sys_oe);
            results_oe = [results_oe; [fit, nb,nf,nk,aic(sys_oe),fpe(sys_oe)]];
            curr_it = curr_it+1;
            waitbar(curr_it/tot_it, f,['OE: ', num2str(round(curr_it/tot_it*100,2)),'%...    Elapsed time: '...
                   ,datestr(seconds(toc),'HH:MM:SS')]);
        end
    end
end
close(f);
save('C:\Users\Joachim\Documents\Prosjekt\Gas lift model\GasLiftModel\Results\results_oe3.mat','results_oe');


%% Finding the best model based on AIC and FPE 
[best_arx, best_arm, best_oe] = find_best_model(results_arx, results_armax, results_oe);

sys_arx = arx(train_data,best_arx(2:4));
figure(1);
compare(test_data,sys_arx); 

sys_arm = armax(train_data, best_arm(2:5));
figure(2);
compare(test_data,sys_arm);

sys_oe = oe(train_data, best_oe(2:4));
figure(3);
compare(test_data,sys_oe);





