%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

function [x_aseb, y_aseb, X_res_Tot, Y_res_Tot] = nodefect(layup_nodefect,layup_defect, Nber_nodefect, Nber_defect, disc_density, t_Hres, x_width, x_al, x_data, y_data)  

% Remove resin thickness on top layer since thickness are defined with thickness
layup_defect(1,1) = layup_defect(1,1) - t_Hres;

%% Total Layup of no defect layers
Tot_layup = layup_nodefect;

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
    
X_data = [l_scarf Xmov_data xend_data];
Y_data = [ybeg_data Ymov_data yend_data]; 
Y_data = linspace(ybeg_data(1,1), yend_data(end), length(Y_data));
X_data = linspace(l_scarf(1,1),  xend_data(end), length(X_data));

% Create matrix
x_Wrinkles = repmat(X_data, Nber_nodefect,1);
y_Wrinkle = repmat(Y_data, Nber_nodefect,1);
y_noWrinkle = y_Wrinkle(:,1).*ones(size(y_Wrinkle));

% Create temporary layer which contains layup information
temp_lays = layup_nodefect.* ones(length(layup_nodefect), length(y_noWrinkle)); 
temp_lay = [zeros(1, length(y_noWrinkle)); temp_lays];
for ia = 1 : Nber_nodefect
    temp_lay(ia + 1,:) = temp_lay(ia + 1,:) +  temp_lay(ia ,:);
    y_Diff(ia,:) = y_noWrinkle(ia,:) - temp_lay(ia+1,:); 
end
y_Wrinkles = flip(y_Diff);   

%% Bottom layup considering info from layup_nodefect
y_asebB    = layup_nodefect.*ones(size(y_Wrinkles));
y_asebB0   = y_asebB(1:end-1,:); % Bottom layer not considered since it will be moved down at y=0
y_tempB    = repmat(y_asebB0(end,:), Nber_nodefect-1,1);
temp_laysB = flip(layup_nodefect(1:end-1,:)).* ones(length(layup_nodefect(1:end-1,:)), length(y_tempB));  
temp_layB  = [ zeros(1, length(y_tempB)); temp_laysB];

for ia = 1 : Nber_nodefect-1
    temp_layB(ia + 1,:) = temp_layB(ia + 1,:) +  temp_layB(ia ,:);
    y_aseb_bot(ia,:)    = y_tempB(ia,:) + temp_layB(ia+1,:) ; 
end  

% Top layup considering info from layup_defect (not taking defect by taking last point from y_aseb_bot and generate layers)
y_Wrinkle = repmat(y_aseb_bot(end,:), Nber_defect,1);
temp_lays = flip(layup_defect).* ones(length(layup_defect), length(y_Wrinkle)); 
temp_lay = [ zeros(1, length(y_Wrinkle)); temp_lays];
for ia = 1 : Nber_defect
    temp_lay(ia + 1,:) = temp_lay(ia + 1,:) +  temp_lay(ia ,:);
    y_aseb_top(ia,:) = y_Wrinkle(ia,:) + temp_lay(ia+1,:) ; 
end

y_aseb = [y_asebB(end,:);y_aseb_bot; y_aseb_top];
y_aseb = y_aseb - y_asebB(end,:);

x_aseb = X_data(1,:).* ones(size(y_aseb)); 

% Discretize according to discretizing density
x_aseb = linspace(x_aseb(1,1), x_aseb(1, end), disc_density);
y_aseb = y_aseb(:,1).* ones(size(x_aseb)); 
x_aseb = x_aseb.* ones(size(y_aseb)); 

X_res_Tot = x_aseb-1; 
Y_res_Tot = y_aseb; 
            
end