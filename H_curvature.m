%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Create curvature using inverse parabolic function

function [Ycurve_aseb] = H_curvature(X_res_Tot, Y_res_Tot, a_curve)

% Take half length of x-axis
xcurveL = X_res_Tot(1,1:length(X_res_Tot(1,:))/2);

% Define a function (parabola in this case) to make curve for half length of x-axis
% Note: higher +a = narrow (flat curve), lower +a= wider (higher curve)
func = sqrt(xcurveL)/sqrt(a_curve); % inverse of a.x^2 parabola 

% Mid point to connect two curves
func_mid = (func(1,end) + func(1,end)*0.001);

% Connect the curves
ycurve_aseb  = [func func_mid flip(func)];

% Add the values Y_res_Tot at ycurve_aseb
Ycurve_aseb = Y_res_Tot + ycurve_aseb(1,1:length(Y_res_Tot));

end