function [best_arx, best_arm, best_oe] = find_best_model(results_arx, results_armax, results_oe)
%% Finding the best model from aic. 
% ARX
[min_arx, id_arx] = min(results_arx(:,5));
[min_arx_fpe, id_arx_fpe] = min(results_arx(:,6));
if id_arx ~= id_arx_fpe
    disp('AIC and FPE for ARX does not give the same best model');
end  
best_arx = results_arx(id_arx,:);

% ARMAX
[min_arm, id_arm] = min(results_armax(:,6));
[min_arm_fpe, id_arm_fpe] = min(results_armax(:,7));
if id_arm ~= id_arm_fpe
    disp('AIC and FPE for ARMAX does not give the same best model');
end  
best_arm = results_armax(id_arm,:);

% OE
[min_oe, id_oe] = min(results_oe(:,5));
[min_oe_fpe, id_oe_fpe] = min(results_oe(:,6));
if id_oe ~= id_oe_fpe
    disp('AIC and FPE for OE does not give the same best model');
end  
best_oe = results_oe(id_oe,:);

end
