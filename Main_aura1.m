%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

clc; fprintf('Beginning to run %s.m.\n', mfilename); close all; clear; 
clear global; workspace; format long g; format compact; 
tic
%% Define material thickness 

% Thickness of resin interface in between layers
t_Hres = 0.08; % [mm]         

% Thickness of UD fibre lamina
t_ud   = 0.8;       ud   = t_ud + t_Hres;    % [mm]  

% Thickness of biax fibre lamina
t_biax = 0.55;      biax = t_biax + t_Hres;  % [mm]  

% Thickness of resin in between repair and parent laminate                                                 
t_Vres = 0.08; % [mm]

% Thicknes of the overlaminate according to material
t_overmat  = biax; % [mm]
     
%% Layup Selection Process                                                                                                     
flag   = 1;                                                                                                                                            
%0 = Composite with defect and layup contains two different lamina materials; 
%1 = Composite with defect and layup contains one type of lamina material; 
%2 = Composite without defect and layup contains two different lamina materials; 
%3 = Composite without defect and layup contains one type of lamina material.

layup_defect   = [ biax biax biax ]';            % [top to bottom]                                                                           
layup_nodefect = [ biax biax biax ]';            % [top to bottom]

%% 0 with same material in repair and parent, 1 with different material
matrep  = 0;

%% Loading: % 0 for tension, 1 for bending
loading   = 0; 
x_disp    = 0.005;  % U1 displacement [mm]
theta_rot = 0.1745; % UR3 rotational displacement [rad]

%% File name of ABAQUS .inp file
namefile = 'Triunfo2';
namefileinp = '.inp';
fileID = fopen(strcat(namefile,namefileinp), 'w'); 

%% Resin Element: 0: Solid Homogenous, 1: Cohesive Elements
resinEl = 1;

%% Overlapping layer length                                                         
xbiax_olap = 50.0;               % [mm]                        
xud_olap   = 100.0;              % [mm]                          
xover_olap = 50.0;               % Overlap length of overlamina [mm]

%% Discretising density (dependent on how defined ABAQUS element)                                                   
disc_density = 1500;                                                                                                                                               

%% Width of defect                  
x_defect = 2;                    % [mm]                                                                                                        
                                                                                                                           
%% Geometry of 2D model                                            
x_width = 500;                   % [mm] width from middle of defect                   
x_al    = 50;                    % [mm] Additional allowance at the edges                 

%% Ply by ply match
mat_change = 0;  % 0: ply by ply match, 1: ply by ply mismatch

%% Overlamina Addition on top of repair
over_lam   = 1;  % 0 no, 1 yes for adding laminate
over_olap  = 1;  % 0 no, 1 yes for adding length on the overlaminate

%% Element number
% NOTE: there is a limit of the combination of number, 
% since it is dependendent on disc_density (if layupdefect is 3, then width
% has 7 inputs, since 3*2 + 1 =7)
width_vec   = [  80 65 45 150 45 65 80 ];    

%% Curve panel
curve = 0;                       % 0 no, 1 yes
a_curve = 200;                   % Parabola curve 
    
%% Guess values for Gaussian parameters                          
amplitude_guess = 8;             % a: also called as global 'c'                    
centre_guess    = 200;           % b: peak position                               
sigma_guess     = 8;             % c: width                                     
                                                                           
%% PART1: Data processing of Defect
layup_nodefect = [layup_nodefect; layup_nodefect(end)];
layup          = [ flip(layup_nodefect); flip(layup_defect) ];      
Nber_defect    = numel(layup_defect);                                                  
Nber_nodefect  = numel(layup_nodefect);                                                                                                                                                                                                                                                  

% Input Data Defect
filename = 'CASE1.xlsx';
sheet    = 1;
Wrinkles = xlsread(filename, sheet);
x_data   = Wrinkles(:,1)'; y_data = Wrinkles(:,2)';
[Wrinkle, Tot_layup, Y_data, X_data] = A_dataprocess(flag, layup_defect,layup_nodefect,Nber_defect, x_width, x_al, x_data, y_data, t_Hres);

%% PART 2: Start loop for defect creation interpolation and create defect according to number of layers with defect
if (flag == 0) || (flag == 1)
[x_store, y_store] = B_peakfitting(flag, Wrinkle, Nber_defect, disc_density, centre_guess, sigma_guess, amplitude_guess, layup_defect, t_Hres);
% Save in order to load later on to decrease running time for high dis_density
% save('Data/Gauss_inputCase1.mat','x_store','y_store')
% load Data/Gauss_inputCase1.mat;

%% PART 3: Create healthy layers below defected layers by extending the layers according to U-def by knowing the last bottom lamina with defect
[x_aseb, y_aseb, x_hty, y_hty] = C_healthylayers(x_store, y_store, layup_nodefect, layup_defect);
    
%% PART 4: Detect defect by locating mid-x_axis of defect, then start defect removal in right section, then take size of removal in right and use it for left section, since it is symmetric  
[xrep_left_aseb, yrep_left_aseb, xrep_right_aseb, yrep_right_aseb,x_right_par, y_right_par, x_left_par, y_left_par, L_remrow] = D_removedefect(layup, flag,x_aseb, y_aseb, x_defect, Nber_nodefect, layup_defect, ud, biax, xud_olap, xbiax_olap);

%% PART 5: Create repair layers by combining left and right 'aseb' matrix
[x_rep, y_rep, Y_top] = E_repairlayers(yrep_left_aseb, xrep_left_aseb, yrep_right_aseb, xrep_right_aseb, x_hty, y_right_par, y_left_par, L_remrow);

%% PART 6: Contruct resin in between parent and repair laminate
[xLsec, xRsec, xRsecEnd, Xad_left, Xad_right] = F_resinParRep(L_remrow, x_left_par, y_left_par,x_right_par, y_right_par, t_Vres, Tot_layup, x_rep, y_rep);

%% PART 7: For element width change to with coarse to fine element widths
Y_res_tot = [y_hty; Y_top];  % assemble y values of all composite
[X_res_Tot, Y_res_Tot] = G_elementchange(width_vec, Tot_layup, xLsec, xRsec, xRsecEnd, Xad_left, Xad_right, Y_res_tot);    

%% PART 8: Adding curvature if requested
if curve == 1 % with curvature in panel
    [Ycurve_aseb] = H_curvature(X_res_Tot, Y_res_Tot, a_curve);

%% PART 9: Tracing all x and y to for parent and resins, etc. It is needed to make nodes and elements        
    [Xres_el, Yres_el, rep_lay, parTrR, parTrL, hty_lay] = I_trace(X_res_Tot, Ycurve_aseb, Xad_left, Xad_right, Nber_defect,Nber_nodefect, flag);

%% PART 10: Adding over laminate
    if (over_lam == 1)
        [x_res_Tot, y_res_Tot, y_over, xy_over] = H_overlaminate(X_res_Tot, Ycurve_aseb, t_overmat, t_Hres, over_olap, rep_lay, xover_olap);
    else if (over_lam == 0)
            x_res_Tot = X_res_Tot; y_res_Tot = Ycurve_aseb;
        end
    end

%% Part 11: Insert interface layer at spanwise horizontal direction
    [X_Tot, Y_Tot, xRes, yRes] = J_resinInterface(t_Hres, x_res_Tot, y_res_Tot, over_lam, t_overmat);

else if curve == 0 % Without curvature
        [Xres_el, Yres_el, rep_lay, parTrR, parTrL, hty_lay] = I_trace(X_res_Tot, Y_res_Tot, Xad_left, Xad_right, Nber_defect,Nber_nodefect, flag);  
    if (over_lam == 1)
        [x_res_Tot, y_res_Tot, y_over, xy_over] = H_overlaminate(X_res_Tot, Y_res_Tot, t_overmat, t_Hres, over_olap, rep_lay, xover_olap);
    else if (over_lam == 0)
            x_res_Tot = X_res_Tot; y_res_Tot = Y_res_Tot;
        end
    end
    % Add Resin interface
    [X_Tot, Y_Tot, xRes, yRes] = J_resinInterface(t_Hres, x_res_Tot, y_res_Tot, over_lam, t_overmat);  
    end
end
  
    else if (flag == 2) || (flag == 3) 
        [x_aseb, y_aseb, X_res_Tot, Y_res_Tot] = nodefect(layup_nodefect,layup_defect, Nber_nodefect, Nber_defect, disc_density, t_Hres, x_width, x_al, x_data, y_data);   
        if curve == 1 % with curvature in panel
            [Ycurve_aseb] = H_curvature(X_res_Tot, Y_res_Tot, a_curve);
            [X_Tot, Y_Tot, xRes, yRes] = J_resinInterface(t_Hres, X_res_Tot, Ycurve_aseb, over_lam, t_overmat);
            else if curve == 0 % Without curvature
                [X_Tot, Y_Tot, xRes, yRes] = J_resinInterface(t_Hres, X_res_Tot, Y_res_Tot, over_lam, t_overmat); 
                end
        end
        end
end

% Figure size
x0     = 50;
y0     = 500;
width  = 2000;
height = 400;

if(flag == 0) || (flag == 2)
    idBiax   = find(layup == biax);
    idUd     = find(layup == ud);
    plotbiax = y_aseb(idBiax,:);
    plotud   = y_aseb(idUd,:);
    figure()
    plot(x_aseb(1, :)/1000, plotbiax/1000, '.r');
    hold on
    plot(x_aseb(1, :)/1000, plotud/1000, '.b');
    xlabel('length / m','Interpreter','latex','Fontsize',13);
    ylabel('thickness / m','Interpreter','latex','Fontsize',13);
    xlim(sort([x_aseb(1)/1000 x_aseb(end)/1000]))
    set(gcf,'position',[x0,y0,width,height])
    hold off
    else if (flag == 1) || (flag == 3)
        if layup(1,1) == biax
            idBiax   = find(layup == biax);
            plotbiax = y_aseb(idBiax,:);
            figure()
            plot(x_aseb(1, :)/1000, plotbiax/1000, '.r');
            legend('Biax')
            xlabel('length / m','Interpreter','latex','Fontsize',13);
            ylabel('thickness / m','Interpreter','latex','Fontsize',13);
            xlim(sort([x_aseb(1)/1000 x_aseb(end)/1000]))
            set(gcf,'position',[x0,y0,width,height])
            else if layup(1,1) == ud
                idUd     = find(layup == ud);
                plotud   = y_aseb(idUd,:);
                figure()
                plot(x_aseb(1, :)/1000, plotud/1000, '.b');
                legend('UD')
                xlabel('length / m','Interpreter','latex','Fontsize',13);
                ylabel('thickness / m','Interpreter','latex','Fontsize',13);
                xlim(sort([x_aseb(1)/1000 x_aseb(end)/1000]))
                set(gcf,'position',[x0,y0,width,height])
                end
        end
    end
end

if (flag==0) || (flag==1)
figure()
plot(x_hty(1, :)/1000, y_hty/1000, '.k'); % healthy layers below
hold on
for i10 = 1:length(L_remrow)
    hold on
    plot(xrep_right_aseb{i10,1}(:,:)/1000, yrep_right_aseb{i10,1}(:,:)/1000, '.r'); % repair layers of left
    hold on
    plot(xrep_left_aseb{i10,1}(:,:)/1000, yrep_left_aseb{i10,1}(:,:)/1000, '.r'); % repair layers of right
    hold on
    plot(x_right_par{i10,1}(:,:)/1000, y_right_par{i10,1}(:,:)/1000, '.k'); % repair layers of left
    hold on
    plot(x_left_par{i10,1}(:,:)/1000, y_left_par{i10,1}(:,:)/1000, '.k'); % repair layers of right
    xlabel('length / m','Interpreter','latex','Fontsize',13);
    ylabel('thickness / m','Interpreter','latex','Fontsize',13);   
    set(gcf,'position',[x0,y0,width,height])
    xlim(sort([x_aseb(1)/1000 x_aseb(end)/1000]))
    hold off
end

figure()
plot(x_hty(1, :)/1000, y_hty/1000, '.k'); % healthy layers below
hold on
for i11 = 1:length(L_remrow)
    hold on
    plot(x_right_par{i11,1}(:,:)/1000, y_right_par{i11,1}(:,:)/1000, '.k'); % parent layers of left
    hold on
    plot(x_left_par{i11,1}(:,:)/1000, y_left_par{i11,1}(:,:)/1000, '.k'); % parent layers of right
    xlabel('length / m','Interpreter','latex','Fontsize',13);
    ylabel('thickness / m','Interpreter','latex','Fontsize',13);
    xlim(sort([x_aseb(1)/1000 x_aseb(end)/1000]))
    set(gcf,'position',[x0,y0,width,height])
    hold off
end

figure()
plot(x_hty(1, :)/1000, y_hty/1000, '.k'); % healthy layers below
hold on
for i11 = 1:length(L_remrow)
    hold on
    plot(x_right_par{i11,1}(:,:)/1000, y_right_par{i11,1}(:,:)/1000, '.k'); % parent layers of left
    hold on
    plot(x_left_par{i11,1}(:,:)/1000, y_left_par{i11,1}(:,:)/1000, '.k'); % parent layers of right
    hold on
    plot(x_rep{i11,1}(:,:)/1000, y_rep{i11,1}(:,:)/1000, 'Color',[0.4660 0.6740 0.1880]); % repair layers
    xlabel('length / m','Interpreter','latex','Fontsize',13);
    ylabel('thickness / m','Interpreter','latex','Fontsize',13);
    xlim(sort([x_aseb(1)/1000 x_aseb(end)/1000]))
    set(gcf,'position',[x0,y0,width,height])
    hold off
end
end

figure()
plot(X_res_Tot(1, :)/1000, yRes/1000, 'Color',[0.9290 0.6940 0.1250], 'LineWidth',3,'MarkerSize',11);
hold on
plot(X_res_Tot(1, :)/1000, Y_res_Tot/1000, '.k', 'LineWidth',3,'MarkerSize',11);
xlabel('length / m','Interpreter','latex','Fontsize',13);
ylabel('thickness / m','Interpreter','latex','Fontsize',13);
xlim(sort([X_res_Tot(1)/1000 X_res_Tot(end)/1000]))
set(gcf,'position',[x0,y0,width,height])

%% Print values
if (flag==0)
    formatSpec = 'Repaired Composite Panel with UD and BIAX material \n';
    fprintf(formatSpec);
else if (flag==1)
        formatSpec = 'Repaired Composite Panel with one material \n';
        fprintf(formatSpec);
    else if (flag==2)
            formatSpec = 'Non-Repaired Composite Panel with UD and BIAX material \n';
            fprintf(formatSpec);
        else if (flag==3)
                formatSpec = 'Non-Repaired Composite Panel with one material \n';
                fprintf(formatSpec);
            end
        end
    end
end

Tot_thick = Y_Tot(end,1);
formatSpec = 'Total thickness of composite is %4.2f mm \n';
fprintf(formatSpec,Tot_thick);

Tot_length = X_Tot(1,end); 
formatSpec = 'Total length of composite is %4.2f mm \n';
fprintf(formatSpec,Tot_length);

if (flag == 0) || (flag == 1) 
    Tot_repair = rep_lay{end,1}(end,1) - rep_lay{end,1}(1,1);
    formatSpec = 'Total length of repair is %4.2f mm \n';
    fprintf(formatSpec,Tot_repair);
end

%% Part 12: ABAQUS Results processing for meshing

% Prepare input for abaqus inp file
LineNberExt = size(X_Tot,1); % y-direction
X_nodenber  = length(X_Tot);  % x-direction
x_extTot    = cell(1, LineNberExt - 1);
y_extTot    = cell(1, LineNberExt - 1);
NberNodef   = Nber_nodefect - 1; % -1 to eliminate the dummy layer
nberrep     = (NberNodef + (NberNodef + 1)); % Number of nodefect layers + Number of resin layers wrt to the next layer
LayerName   = (0:LineNberExt-1)';
ResinName   = (0:2:LineNberExt-1)';
LayupMat    = [flip(layup_nodefect(1:end-1,:)); flip(layup_defect)];
countLay    = 1;

for i11 = 1 : LineNberExt
    x_extTot{1, countLay} = X_Tot(i11, :); % bottom to top
    y_extTot{1, countLay} = Y_Tot(i11, :); % bottom to top
    countLay = countLay + 1;
end

% Nodes
x_mat = cell2mat(x_extTot)'/1000; % convert cell data to matrix
y_mat = cell2mat(y_extTot)'/1000;
TotNode  = (LineNberExt) * X_nodenber; % Total node number
NodeNber = cell(TotNode, 1);

% Text formatting
FormatY = 2; FormatX = 4;  % 2 digits for layer number , 4 digits for number of points depending on disc_density

% Syntax
NberFormaty = ['%0', num2str(FormatY), '.0f']; NberFormatx = ['%0', num2str(FormatX), '.0f'];

% Assign node number to each node 
countNber = 1;
for i6 = 1:LineNberExt
    for i7 = 1 : X_nodenber
        NodeNber{countNber, 1} = [num2str(i6, NberFormaty), num2str(i7, NberFormatx)];
        countNber = countNber + 1;
    end
end

%% Print Abaqus input file according to case

% With same material in repair and parent laminate
if (matrep == 0)
    % With repair case
    if (flag == 0) || (flag == 1)
        save inputwithrepair.mat
            if (loading == 0) % Tension Case
                if (resinEl == 1)
                    run('K_abaqusrepairCZM.m'); % Resin modelled with cohesive elements
                else if (resinEl == 0)
                        run('K_abaqusrepairSolid.m'); % Resin modelled with solid elements
                    end
                end
            else if (loading == 1) % Bending case
                if (resinEl == 1)
                    run('L_abaqusrepairCZM.m'); % Resin modelled with cohesive elements 
                else if (resinEl == 0)
                        run('L_abaqusrepairSolid.m'); % Resin modelled with solid elements
                    end
                end 
                end
            end
    end
    
    % Without repair case
    if (flag == 2) || (flag == 3)
        save inputwithoutrepair.mat
            if (loading == 0) % Tension Case
                if  (resinEl == 1)
                    run('K_abaquspristineCZM.m'); % Resin modelled with cohesive elements
                else if (resinEl == 0)
                        run('K_abaquspristineSolid.m'); % Resin modelled with solid elements
                        end
                end
                else if (loading == 1) % Bending case 
                        if  (resinEl == 1)
                            run('L_abaquspristineCZM.m'); % Resin modelled with cohesive elements
                        else if (resinEl == 0)
                                run('L_abaquspristineSolid.m'); % Resin modelled with solid elements
                            end
                        end
                    end
            end
    end
    
    % With Different material in repair
    else if (matrep == 1)
            save inputwithrepair.mat
            if (loading == 0) % Tension Case
                if (resinEl == 1)
                    run('K_abaqusmatdiffCZM.m'); % Resin modelled with cohesive elements
                else if (resinEl == 0)
                        run('K_abaqusmatdiffSolid.m'); % Resin modelled with solid elements
                    end
                end
            else if (loading == 1) % Bending case
                if (resinEl == 1)
                    run('L_abaqusmatdiffCZM.m'); % Resin modelled with cohesive elements 
                else if (resinEl == 0)
                        run('L_abaqusmatdiffSolid.m'); % Resin modelled with solid elements
                    end
                end 
                end
            end
        end
end
%% Until here the code can be run with MATLAB, after this moment a full ABAQUS license it's requiered. 
%% Write the necessary lines in .inp so that ABAQUS gives the dessired ouput
% Change the file so that microsoft nopad can read the lines correctly \r\n
% instead of \n and add the necessary lines for ABAQUS
% Delete the overlapping lines for the cohesive file
% if resinEl==1
%     func_replace_string2(strcat(namefile,namefileinp), strcat(namefile,namefileinp), '', '');
% end
func_replace_string2(strcat(namefile,namefileinp), strcat(namefile,namefileinp), '*End Step', '*NODE PRINT, nset=RP_LEFT, frequency=1')
fid = fopen(strcat(namefile,namefileinp),'a');
if loading == 0
    fprintf(fid,'RF1,RF2\r\n');
elseif loading == 1
    fprintf(fid,'RM3,RF2\r\n');
end
fprintf(fid, '*End Step\r\n'); 
fclose(fid);
%% Send the .inp file to ABAQUS, submit and read the output
if loading == 0
    [RF1,RF2] = Run_job_request_outputs_function(namefile); 
elseif loading == 1
    [RM3,RF2] = Run_job_request_outputs_function(namefile);
end
%% Calculate the Young's modulus or the Bending stiffness
if loading == 0
    E = Youngs_modulus(x_disp,RF1(10),Tot_length,Tot_thick);
elseif loading == 1
    D = Bending_stiffness(theta_rot,RM3(10),Tot_length);
end
toc





