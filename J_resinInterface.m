%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: Add resin interface between top and bottom layers

function [X_Tot, Y_Tot, xRes, yRes] = J_resinInterface(t_Hres, X_res_Tot, Y_res_Tot, over_lam, t_overmat)
       
LNberExt = size(X_res_Tot, 1);
LineNberRes = LNberExt - 1;
LineNberTot = LNberExt + LineNberRes;

for i13 = 1:LineNberRes-1
    % Make adhesive layers line, NOTE all layers already contains adhesive horizontaly, this making elements for them
    yRes(i13,:) = Y_res_Tot(end-i13,:) - t_Hres;
    xRes(i13,:) = X_res_Tot(i13,:);
end

% Assemble all layers with ahesive
if over_lam == 1
    Y_over = Y_res_Tot(end,:) - t_overmat  + t_Hres;
    Y_Tot = sortrows([yRes; Y_res_Tot]);
    X_Tot = [xRes; X_res_Tot];
else
    Y_Tot = sortrows([yRes; Y_res_Tot]);
    X_Tot = [xRes; X_res_Tot];
end


end