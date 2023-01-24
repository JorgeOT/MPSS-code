%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

%Description: Take input file which contains information of one layer
%with wrinkle. This defect will be adjusted according to how much the
%length of the scarfing will be performed during repair, and also how many
%layers are with defect.

function [Wrinkle, Tot_layup, Y_data, X_data] = A_dataprocess(flag, layup_defect,layup_nodefect,Nber_defect, x_width, x_al, x_data, y_data, t_Hres)

% Compute total thickness of layup
if (flag == 0) 
    layup_defect(1,1) = layup_defect(1,1) - t_Hres;
    Tot_layup = [layup_defect; layup_nodefect];
else if (flag == 1)
        Tot_layup = [layup_defect; layup_nodefect];
        else if (flag == 2) || (flag == 3)
            Tot_layup = layup_nodefect;                                             
            end
    end
end

% Moving data
x_scarf = x_width + x_al;                                                   % Adding allowance at the end
y_scarf = sum(Tot_layup);                                                   % [mm] Total thickness of composite , adhesive in bottom is included

Xmov_data = x_data + x_scarf;
Ymov_data = y_data + y_scarf;

% Preparing for rigth and left ends of the defect
l_scarf = linspace(1, x_scarf, x_scarf);
xend_data = Xmov_data(1,end) + l_scarf;
yend_data = Ymov_data(1,end)*ones(size(l_scarf));
ybeg_data = Ymov_data(1,1)*ones(size(l_scarf));

% Assemble manipulated data
X_data = [l_scarf Xmov_data xend_data];
Y_data = [ybeg_data Ymov_data yend_data];

% Storing into cell
Wrinkle{1,1}(:,1) = X_data;
Wrinkle{1,1}(:,2) = Y_data;
Wrinkle = flip(Wrinkle,2);

end