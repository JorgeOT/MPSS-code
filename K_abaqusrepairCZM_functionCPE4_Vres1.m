%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)
function [] = K_abaqusrepairCZM_functionCPE4_Vres1(V12RESIN,E1RESIN,I_o,I_R,Damstab,t_Vres,flag, layup,t_Hres,t_overmat,E1UD,E1BIAX,E2UD,E2BIAX,E3UD,E3BIAX,V12UD,V12BIAX,V13UD,V13BIAX,V23UD,V23BIAX,G12UD,G12BIAX,G13UD,G13BIAX,G23UD,G23BIAX,first_Step,Tperiod,Min_incr,Max_incr,Viscosity,Daminit_UD,Daminit_UD_DIFF,Daminit_BIAX,Daminit_BIAX_DIFF,Daminit_Resin1,Daminit_Resin2,Daminit_Resin3,Damev_UD,Damev_BIAX,Damev_Resin,Frac_UD1,Frac_UD2,Frac_UD3,Frac_BIAX1,Frac_BIAX2,Frac_BIAX3,Frac_Resin1,Frac_Resin2,Frac_Resin3,Frac_UD1_DIFF,Frac_UD2_DIFF,Frac_UD3_DIFF,Frac_BIAX1_DIFF,Frac_BIAX2_DIFF,Frac_BIAX3_DIFF,Knn,Kss,Ktt)
% Load workspace from main
load('inputwithrepair.mat');

%% Nodes and Part
LineText = '** ------------------------------------';
fprintf(fileID, '%s\r\n', LineText);
LineText = '**';
fprintf(fileID, '%s\r\n', LineText);
LineText = '** PARTS';
fprintf(fileID, '%s\r\n', LineText);
LineText = '**';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Part, name=RepairedPanel';
fprintf(fileID, '%s\r\n', LineText);
LineText = '** NodeNber	XCoord	YCoord';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText);

for i8 = 1:TotNode
	LineText = [NodeNber{i8,1},',', num2str(x_mat(i8, 1), '%02.8f'),',',num2str(y_mat(i8, 1), '%02.8f') ];
    fprintf(fileID, '%s\r\n', LineText);
end

%% Elements
% Quadratic element and dimension
d2 = 4; % 2d element, element type, plane strain 'CPE4', plane stress 'CPS4' 
TotElt = (LineNberExt - 1)*(X_nodenber - 1); % Total element number 
EltName = cell(TotElt, 1);
EltNodes = cell(TotElt, d2); 
CompId = zeros(TotElt, 1); 
countElt = 1;

% Assign element number and Node order
for i9 = 1 : LineNberExt - 1
    for i10 = 1 : X_nodenber - 1
        % Composite or not? 
        if mod(i9, 2) ~= 0 % odd
           CompId(countElt, 1) = 1; % Composite = 1, Adhesive = 0   
        end        
        % Element name y-direction
        EltName{countElt, 1} = [num2str(i9, NberFormaty), num2str(i10, NberFormatx)];
        
        % Node 1
        EltNodes{countElt, 1} = [num2str(i9, NberFormaty), num2str(i10, NberFormatx)];
        % Node 2
        EltNodes{countElt, 2} = [num2str(i9, NberFormaty), num2str(i10 + 1, NberFormatx)];
        % Node 3
        EltNodes{countElt, 3} = [num2str(i9 + 1, NberFormaty), num2str(i10 + 1, NberFormatx)];
        % Node 4
        EltNodes{countElt, 4} = [num2str(i9 + 1, NberFormaty), num2str(i10, NberFormatx)];
        countElt = countElt + 1;
    end
end

%% Sets: Resin Interface
% Combine x_mat and y_mat, and locate the resin elemeent (to be combined)
mat_com = [ x_mat   y_mat];
dif_el = str2num([num2str(1, NberFormaty), num2str(0, NberFormatx)]);

if over_olap == 1
    TotEltres = (LineNberExt - 3)*(X_nodenber - 1); % Total element number to not consider last interface layer
else
    TotEltres = (LineNberExt - 1)*(X_nodenber - 1); % Total element number 
end

LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Element, type=COH2D4, Elset=Einterface'; % Cohesive elements
fprintf(fileID, '%s\r\n', LineText); 
for ii11 = 1 : TotEltres
    if CompId(ii11, 1) == 0
        LineText = EltName{ii11, 1};
            for ii12 = 1 : size(EltNodes, 2)
                LineText = [LineText,', ',	EltNodes{ii11, ii12}];
            end
        fprintf(fileID, '%s\r\n', LineText);
    end
end

% Add elements and nodes of the overlaminate with varied length
if over_olap == 1
    loc_overres = ismember(mat_com, xy_over/1000, 'rows');
    Node_overres = NodeNber(loc_overres);
    for m = 1:length(Node_overres)-1
        Elt_Over = [num2str(str2num(Node_overres{m,:}) - (dif_el*2)), ', ', num2str(str2num(Node_overres{m,:}) - (dif_el*2)), ', ' ,num2str(str2num(Node_overres{m,:}) - (dif_el*2) + 1) , ', ',  num2str(str2num(Node_overres{m,:})- dif_el + 1), ', ', num2str(str2num(Node_overres{m,:}) - dif_el) ];
        fprintf(fileID, '%s\r\n', Elt_Over);
    end
end

%% Set: Overlaminate material   
if over_lam == 1
    loc_over = ismember(mat_com, xy_over/1000, 'rows');
    Node_over = NodeNber(loc_over);

    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=Eoverlam';
    fprintf(fileID, '%s\r\n', LineTextloop);
    for m = 1:length(Node_over) - 1
        Elt_Over = [num2str(str2num(Node_over{m,:}) - dif_el), ', ', num2str(str2num(Node_over{m,:}) - dif_el), ', ' ,num2str(str2num(Node_over{m,:}) - dif_el + 1) , ', ',  num2str(str2num(Node_over{m,:}) + 1), ', ', num2str(str2num(Node_over{m,:})) ];
        fprintf(fileID, '%s\r\n', Elt_Over);
    end
end

%% Set: Resin in between parent and repair
% Combine x_mat and y_mat, and locate the resin elemeent (to be combined)
res_com = [ Xres_el/1000 Yres_el/1000 ];

% Locate res_com at mat_com and get indeces to get the Nodenumber
locres_com = ismembertol(mat_com,res_com,'ByRows',true);

% Take corresponfing nodes
Node_Resloc = NodeNber(locres_com);

% Corresponding 4th and 3rd node
Node_Res = [Node_Resloc(1:2:end) Node_Resloc(2:2:end)];  

% Divide by left and right values
Node_ResL = Node_Res(1:2:end,:);
Node_ResR = Node_Res(2:2:end,:);

% Take idx layers of odd composites layers
% comp     = (nberrep+2 : 2 : LineNberExt)'; % all lamina layers has name in odd numbers
comp     = (nberrep : 2 : LineNberExt)'; % all lamina layers has name in odd numbers
comp_mat = flip(layup_defect);

% Build elements (1st and 2nd) and combine elem pairs.
for ii8 = 1:length(Node_ResL)
    % Make element name formatx
    repResmat_NameA = num2str(comp(ii8,1));

    % left resin
    repResLmat_Name1 = Node_ResL{ii8,1}(FormatY+1:end); % 4th
    repResLmat_Name2 = Node_ResL{ii8,2}(FormatY+1:end); % 3rd
    repResLmat_1{ii8, 1} = str2num(strcat(repResmat_NameA,repResLmat_Name1)); % 1st node
    repResLmat_2{ii8, 1} = str2num(strcat(repResmat_NameA,repResLmat_Name2)); % 2nd node
    repResLmat(ii8,:)   = [repResLmat_1{ii8, 1} repResLmat_2{ii8, 1}];

    % right resin
    repResRmat_Name1 = Node_ResR{ii8,1}(FormatY+1:end); % 4th
    repResRmat_Name2 = Node_ResR{ii8,2}(FormatY+1:end); % 3rd
    repResRmat_1{ii8, 1} = str2num(strcat(repResmat_NameA,repResRmat_Name1)); % 1st node
    repResRmat_2{ii8, 1} = str2num(strcat(repResmat_NameA,repResRmat_Name2)); % 2nd node
    repResRmat(ii8,:)   = [repResRmat_1{ii8, 1} repResRmat_2{ii8, 1}];
end

% Combine in one matrix
repResmat = [repResLmat; repResRmat];

% Print
LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Element, type=COH2D4, Elset=Eresin';
fprintf(fileID, '%s\r\n', LineText); 
for ii8 = 1:length(Node_Res)
    % Create element
    Elt_Res = [num2str(repResmat(ii8,1)), ', ', num2str(repResmat(ii8,1)), ', ' ,num2str(repResmat(ii8,2)) , ', ',  num2str(repResmat(ii8,2) + dif_el), ', ', num2str(repResmat(ii8,1) + dif_el) ];
    fprintf(fileID, '%s\r\n', Elt_Res);
end

%% Sets: Repair layers
% Locate res_com at mat_com and get indeces to get the Nodenumber

if mat_change == 1  
    layup_change = [layup_defect(2:end); t_Hres];
    
    % Resin richness
    locrep_res = ismember(mat_com,rep_lay{1,1}/1000, 'rows');
    Node_Represloc = NodeNber(locrep_res);
    
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=Eresinrich';
    fprintf(fileID, '%s\r\n', LineTextloop);
    for m = 1:length(Node_Represloc)
        Elt_Over = [num2str(str2num(Node_Represloc{m,:}) - (dif_el*2)), ', ', num2str(str2num(Node_Represloc{m,:})- (dif_el*2)), ', ' ,num2str(str2num(Node_Represloc{m,:}) - (dif_el*2) + 1) , ', ',  num2str(str2num(Node_Represloc{m,:}) - (dif_el) + 1), ', ', num2str(str2num(Node_Represloc{m,:})- (dif_el)) ];
        fprintf(fileID, '%s\r\n', Elt_Over);
    end
    
    % Composite repair layers
    for i = 1:length(layup_change)-1
        if layup_change(i,:) == ud
            loc_changeUd(:,1) = ismember(mat_com,rep_lay{i+1,1}/1000, 'rows');
            Node_changeUd{i,1} = NodeNber(loc_changeUd(:,1)); 
            else if layup_change(i,:) == biax
                loc_changeBiax(:,1) = ismember(mat_com,rep_lay{i+1,1}/1000, 'rows');
                Node_changeBIAX{i,1} = NodeNber(loc_changeBiax(:,1));
            end
        end
    end

    if any(layup_change == ud)
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EREPAIRUD1';%ErepairUD1
    fprintf(fileID, '%s\r\n', LineTextloop);
    for i = 1:length(Node_changeUd)
        if isempty(Node_changeUd{i,1}) == 0
        for m = 1:length(Node_changeUd{i,1})
            repRepUD_NameA      = num2str(comp(i+1,1));
            repRepUD_NameB_temp = Node_changeUd{i,1}; % 1st
            repRepUD_NameB = repRepUD_NameB_temp{m,1}(FormatY+1:end);
            Node_changeUdloc{i, m} = str2num(strcat(repRepUD_NameA,repRepUD_NameB)); % 1st node, Element name
            Elt_Over = [num2str(Node_changeUdloc{i, m}), ', ', num2str(Node_changeUdloc{i, m}), ', ' ,num2str(Node_changeUdloc{i, m}+1) , ', ',  num2str(Node_changeUdloc{i, m} +1 + dif_el), ', ', num2str(Node_changeUdloc{i, m} + dif_el) ];
            fprintf(fileID, '%s\r\n', Elt_Over);
        end
        else
            continue;
        end
    end
    end
    
    if any(layup_change == biax)    
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EREPAIRBIAX1';%ErepairBIAX1
    fprintf(fileID, '%s\r\n', LineTextloop);
    for i = 1:length(Node_changeBIAX)
        if isempty(Node_changeBIAX{i,1}) == 0
            for m = 1:length(Node_changeBIAX{i,1})
            repRepBIAX_NameA      = num2str(comp(i+1,1));
            repRepBIAX_NameB_temp = Node_changeBIAX{i,1}; % 1st
            repRepBIAX_NameB = repRepBIAX_NameB_temp{m,1}(FormatY+1:end);
            Node_changeBIAXloc{i, m} = str2num(strcat(repRepBIAX_NameA,repRepBIAX_NameB)); % 1st node, Element name
            Elt_Over = [num2str(Node_changeBIAXloc{i, m}), ', ', num2str(Node_changeBIAXloc{i, m}), ', ' ,num2str(Node_changeBIAXloc{i, m}+1) , ', ',  num2str(Node_changeBIAXloc{i, m} +1 + dif_el), ', ', num2str(Node_changeBIAXloc{i, m} + dif_el) ];
            fprintf(fileID, '%s\r\n', Elt_Over);
        end
        else
            continue;
        end
    end    
    end
    
else   
for i = 1:Nber_defect
    locrep_com(:,1) = ismember(mat_com,rep_lay{i,1}/1000, 'rows');
    Node_Reploc{i,1} = NodeNber(locrep_com(:,1));
end
    
if (flag==0)
    % Divide by material the repair layers
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EREPAIRUD';%ErepairUD
    fprintf(fileID, '%s\r\n', LineTextloop); 
    for m = 1:length(Node_Reploc)
        for n  = 1:length(Node_Reploc{m,1})
            % Layer name
            if comp_mat(m,1) == ud
                repRepUD_NameA      = num2str(comp(m,1));
                repRepUD_NameB_temp = Node_Reploc{m,1}; % 1st
                repRepUD_NameB = repRepUD_NameB_temp{n,1}(FormatY+1:end);
                repRepUD{n, m} = str2num(strcat(repRepUD_NameA,repRepUD_NameB)); % 1st node, Element name
                Elt_RepUD = [num2str(repRepUD{n, m}), ', ', num2str(repRepUD{n, m}), ', ' ,num2str(repRepUD{n, m}+1) , ', ',  num2str(repRepUD{n, m} +1 + dif_el), ', ', num2str(repRepUD{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_RepUD);
            end
        end
    end
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EREPAIRBIAX';%ErepairBiax
    fprintf(fileID, '%s\r\n', LineTextloop); 
    for m = 1:length(Node_Reploc)
        for n  = 1:length(Node_Reploc{m,1})
            % Layer name
            if comp_mat(m,1) == biax
                repRepBiax_NameA      = num2str(comp(m,1));
                repRepBiax_NameB_temp = Node_Reploc{m,1}; % 1st
                repRepBiax_NameB = repRepBiax_NameB_temp{n,1}(FormatY+1:end);
                repRepBiax{n, m} = str2num(strcat(repRepBiax_NameA,repRepBiax_NameB)); % 1st node, Element name
                Elt_RepBiax = [num2str(repRepBiax{n, m}), ', ', num2str(repRepBiax{n, m}), ', ' ,num2str(repRepBiax{n, m}+1) , ', ',  num2str(repRepBiax{n, m} +1 + dif_el), ', ', num2str(repRepBiax{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_RepBiax);
            end
        end
    end   
else if (flag==1) && layup(1,1) == ud
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EREPAIRUD';%ErepairUD
    fprintf(fileID, '%s\r\n', LineTextloop); 
    for m = 1:length(Node_Reploc)
        for n  = 1:length(Node_Reploc{m,1})
            % Layer name
            if comp_mat(m,1) == ud
                repRepUD_NameA      = num2str(comp(m,1));
                repRepUD_NameB_temp = Node_Reploc{m,1}; % 1st
                repRepUD_NameB = repRepUD_NameB_temp{n,1}(FormatY+1:end);
                repRepUD{n, m} = str2num(strcat(repRepUD_NameA,repRepUD_NameB)); % 1st node, Element name
                Elt_RepUD = [num2str(repRepUD{n, m}), ', ', num2str(repRepUD{n, m}), ', ' ,num2str(repRepUD{n, m}+1) , ', ',  num2str(repRepUD{n, m} +1 + dif_el), ', ', num2str(repRepUD{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_RepUD);
            end
        end
    end    
    else if (flag==1) && layup(1,1) == biax
        LineText = '**--------------------------------------'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineTextloop = '*Element, type=CPE4, Elset=EREPAIRBIAX';%ErepairUD
        fprintf(fileID, '%s\r\n', LineTextloop); 
        for m = 1:length(Node_Reploc)
            for n  = 1:length(Node_Reploc{m,1})
                % Layer name
                if CompId(m,1) == 1
                    repRepBiax_NameA      = num2str(comp(m,1));
                    repRepBiax_NameB_temp = Node_Reploc{m,1}; % 1st
                    repRepBiax_NameB = repRepBiax_NameB_temp{n,1}(FormatY+1:end);
                    repRepBiax{n, m} = str2num(strcat(repRepBiax_NameA,repRepBiax_NameB)); % 1st node, Element name
                    Elt_RepBiax = [num2str(repRepBiax{n, m}), ', ', num2str(repRepBiax{n, m}), ', ' ,num2str(repRepBiax{n, m}+1) , ', ',  num2str(repRepBiax{n, m} +1 + dif_el), ', ', num2str(repRepBiax{n, m} + dif_el) ];
                    fprintf(fileID, '%s\r\n', Elt_RepBiax);
                end
            end
        end
        end
    end
end
end

%% Sets: Parent layers and healthy layers

% Locate parent at mat_com and get indeces to get the Nodenumber
for i = 1:Nber_defect
    locpar_TrR(:,1) = ismember(mat_com, parTrR{i,1}/1000, 'rows');
    locpar_TrL(:,1) = ismember(mat_com, parTrL{i,1}/1000, 'rows');
    Node_parTrR{i,1} = NodeNber(locpar_TrR(:,1));
    Node_parTrL{i,1} = NodeNber(locpar_TrL(:,1));
end

% Locate healthy layers
comp_hty = flip(layup_nodefect(1:end-1,:));

for i = 1:Nber_nodefect-1
    loc_hty(:,1)  = ismember(mat_com, hty_lay{1,i}/1000, 'rows');
    Node_hty{i,1} = NodeNber(loc_hty(:,1));
end

if (flag==0)
    % UD Parent
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EparentUD';
    fprintf(fileID, '%s\r\n', LineTextloop); 
    for m = 1:length(Node_parTrL)
        for n  = 1:length(Node_parTrL{m,1})
            % Layer name
            if comp_mat(m,1) == ud
                repParLUd_NameA = num2str(comp(m,1));
                repParLUd_NameB_temp = Node_parTrL{m,1}; % 1st
                repParLUd_NameB = repParLUd_NameB_temp{n,1}(FormatY+1:end);
                repParLUd{n, m} = str2num(strcat(repParLUd_NameA,repParLUd_NameB)); % 1st node, Element name
                Elt_ParLUd = [num2str(repParLUd{n, m}), ', ', num2str(repParLUd{n, m}), ', ' ,num2str(repParLUd{n, m}+1) , ', ',  num2str(repParLUd{n, m} +1 + dif_el), ', ', num2str(repParLUd{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_ParLUd);
            end
        end
    end
    for m = 1:length(Node_parTrR)
        for n  = 1:length(Node_parTrR{m,1})
            % Layer name
            if comp_mat(m,1) == ud
                repParRUd_NameA = num2str(comp(m,1));
                repParRUd_NameB_temp = Node_parTrR{m,1}; % 1st
                repParRUd_NameB = repParRUd_NameB_temp{n,1}(FormatY+1:end);        
                % Ensure end value is part of nodenumber, if not, it is omitted
                if repParRUd_NameB == NodeNber{end,1}(FormatY+1:end)
                     n = n-1;
                     continue;
                end
                repParRUd{n, m} = str2num(strcat(repParRUd_NameA,repParRUd_NameB)); % 1st node, Element name
                Elt_ParRUd = [num2str(repParRUd{n, m}), ', ', num2str(repParRUd{n, m}), ', ' ,num2str(repParRUd{n, m}+1) , ', ',  num2str(repParRUd{n, m} +1 + dif_el), ', ', num2str(repParRUd{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_ParRUd);
            end
        end
    end
    for m = 1:length(Node_hty)
        for n  = 1:length(Node_hty{m,1})
            % Layer name
            if comp_hty(m,1) == ud
                htyUd = cell2mat(Node_hty{m, 1}(n,1));
                htyUd_Name = htyUd(FormatY+1:end);                
                if htyUd_Name == NodeNber{end,1}(FormatY+1:end)
                     n = n-1; % skip the specific iteration
                     continue;
                end
                Elt_hty = [cell2mat(Node_hty{m, 1}(n,1)), ', ', cell2mat(Node_hty{m, 1}(n,1)), ', ' ,num2str(str2double(Node_hty{m,1}(n,1)) +1) , ', ',  num2str(str2double(Node_hty{m, 1}(n,1)) +1 + dif_el), ', ', num2str(str2double(Node_hty{m, 1}(n,1)) + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_hty);
            end
        end
    end
    % Biax Parent
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineTextloop = '*Element, type=CPE4, Elset=EparentBiax';
    fprintf(fileID, '%s\r\n', LineTextloop); 
    for m = 1:length(Node_parTrL)
        for n  = 1:length(Node_parTrL{m,1})
            % Layer name
            if comp_mat(m,1) == biax
                repParLBiax_NameA = num2str(comp(m,1));
                repParLBiax_NameB_temp = Node_parTrL{m,1}; % 1st
                repParLBiax_NameB = repParLBiax_NameB_temp{n,1}(FormatY+1:end);
                repParLBiax{n, m} = str2num(strcat(repParLBiax_NameA,repParLBiax_NameB)); % 1st node, Element name
                Elt_ParLBiax = [num2str(repParLBiax{n, m}), ', ', num2str(repParLBiax{n, m}), ', ' ,num2str(repParLBiax{n, m}+1) , ', ',  num2str(repParLBiax{n, m} +1 + dif_el), ', ', num2str(repParLBiax{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_ParLBiax);
            end
        end
    end
    for m = 1:length(Node_parTrR)
        for n  = 1:length(Node_parTrR{m,1})
            % Layer name
            if comp_mat(m,1) == biax
                repParRBiax_NameA = num2str(comp(m,1));
                repParRBiax_NameB_temp = Node_parTrR{m,1}; % 1st
                repParRBiax_NameB = repParRBiax_NameB_temp{n,1}(FormatY+1:end);          
                % Ensure end value is part of nodenumber, if not, it is omitted
                if repParRBiax_NameB == NodeNber{end,1}(FormatY+1:end)
                     n = n-1; % skip the specific iteration
                     continue;
                end
                repParRBiax{n, m} = str2num(strcat(repParRBiax_NameA,repParRBiax_NameB)); % 1st node, Element name
                Elt_ParRBiax = [num2str(repParRBiax{n, m}), ', ', num2str(repParRBiax{n, m}), ', ' ,num2str(repParRBiax{n, m}+1) , ', ',  num2str(repParRBiax{n, m} +1 + dif_el), ', ', num2str(repParRBiax{n, m} + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_ParRBiax);
            end
        end
    end
    for m = 1:length(Node_hty)
        for n  = 1:length(Node_hty{m,1})
            % Layer name
            if comp_hty(m,1) == biax
                htyBiax = cell2mat(Node_hty{m, 1}(n,1));
                htyBiax_Name = htyBiax(FormatY+1:end);
                if htyBiax_Name == NodeNber{end,1}(FormatY+1:end)
                     n = n-1; % skip the specific iteration
                     continue;
                end
                Elt_hty = [cell2mat(Node_hty{m, 1}(n,1)), ', ', cell2mat(Node_hty{m, 1}(n,1)), ', ' ,num2str(str2double(Node_hty{m,1}(n,1)) +1) , ', ',  num2str(str2double(Node_hty{m, 1}(n,1)) +1 + dif_el), ', ', num2str(str2double(Node_hty{m, 1}(n,1)) + dif_el) ];
                fprintf(fileID, '%s\r\n', Elt_hty);
            end
        end
    end
    
% FLAG 1 FOR UD
else if (flag==1) && layup(1,1) == ud   
        % UD Parent
        LineText = '**--------------------------------------'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineTextloop = '*Element, type=CPE4, Elset=EparentUD';
        fprintf(fileID, '%s\r\n', LineTextloop); 
        for m = 1:length(Node_parTrL)
            for n  = 1:length(Node_parTrL{m,1})
                % Layer name
                if comp_mat(m,1) == ud
                    repParLUd_NameA = num2str(comp(m,1));
                    repParLUd_NameB_temp = Node_parTrL{m,1}; % 1st
                    repParLUd_NameB = repParLUd_NameB_temp{n,1}(FormatY+1:end);
                    repParLUd{n, m} = str2num(strcat(repParLUd_NameA,repParLUd_NameB)); % 1st node, Element name
                    Elt_ParLUd = [num2str(repParLUd{n, m}), ', ', num2str(repParLUd{n, m}), ', ' ,num2str(repParLUd{n, m}+1) , ', ',  num2str(repParLUd{n, m} +1 + dif_el), ', ', num2str(repParLUd{n, m} + dif_el) ];
                    fprintf(fileID, '%s\r\n', Elt_ParLUd);
                end
            end
        end
        for m = 1:length(Node_parTrR)
            for n  = 1:length(Node_parTrR{m,1})
                % Layer name
                if comp_mat(m,1) == ud
                    repParRUd_NameA = num2str(comp(m,1));
                    repParRUd_NameB_temp = Node_parTrR{m,1}; % 1st
                    repParRUd_NameB = repParRUd_NameB_temp{n,1}(FormatY+1:end);        
                    % Ensure end value is part of nodenumber, if not, it is omitted
                    if repParRUd_NameB == NodeNber{end,1}(FormatY+1:end)
                         n = n-1;
                         continue;
                    end
                    repParRUd{n, m} = str2num(strcat(repParRUd_NameA,repParRUd_NameB)); % 1st node, Element name
                    Elt_ParRUd = [num2str(repParRUd{n, m}), ', ', num2str(repParRUd{n, m}), ', ' ,num2str(repParRUd{n, m}+1) , ', ',  num2str(repParRUd{n, m} +1 + dif_el), ', ', num2str(repParRUd{n, m} + dif_el) ];
                    fprintf(fileID, '%s\r\n', Elt_ParRUd);
                end
            end
        end
        for m = 1:length(Node_hty)
            for n  = 1:length(Node_hty{m,1})
                % Layer name
                if comp_hty(m,1) == ud
                    htyUd = cell2mat(Node_hty{m, 1}(n,1));
                    htyUd_Name = htyUd(FormatY+1:end);
                    if htyUd_Name == NodeNber{end,1}(FormatY+1:end)
                         n = n-1; % skip the specific iteration
                         continue;
                    end
                    Elt_hty = [cell2mat(Node_hty{m, 1}(n,1)), ', ', cell2mat(Node_hty{m, 1}(n,1)), ', ' ,num2str(str2double(Node_hty{m,1}(n,1)) +1) , ', ',  num2str(str2double(Node_hty{m, 1}(n,1)) +1 + dif_el), ', ', num2str(str2double(Node_hty{m, 1}(n,1)) + dif_el) ];
                    fprintf(fileID, '%s\r\n', Elt_hty);
                end
            end
        end
        
% FLAG 1 FOR BIAX
        else if (flag==1) && layup(1,1) == biax
            % Biax Parent
            LineText = '**--------------------------------------'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineTextloop = '*Element, type=CPE4, Elset=EparentBiax';
            fprintf(fileID, '%s\r\n', LineTextloop); 
            for m = 1:length(Node_parTrL)
                for n  = 1:length(Node_parTrL{m,1})
                    % Layer name
                    if comp_mat(m,1) == biax
                        repParLBiax_NameA = num2str(comp(m,1));
                        repParLBiax_NameB_temp = Node_parTrL{m,1}; % 1st
                        repParLBiax_NameB = repParLBiax_NameB_temp{n,1}(FormatY+1:end);
                        repParLBiax{n, m} = str2num(strcat(repParLBiax_NameA,repParLBiax_NameB)); % 1st node, Element name
                        Elt_ParLBiax = [num2str(repParLBiax{n, m}), ', ', num2str(repParLBiax{n, m}), ', ' ,num2str(repParLBiax{n, m}+1) , ', ',  num2str(repParLBiax{n, m} +1 + dif_el), ', ', num2str(repParLBiax{n, m} + dif_el) ];
                        fprintf(fileID, '%s\r\n', Elt_ParLBiax);
                    end
                end
            end
            for m = 1:length(Node_parTrR)
                for n  = 1:length(Node_parTrR{m,1})
                    % Layer name
                    if comp_mat(m,1) == biax
                        repParRBiax_NameA = num2str(comp(m,1));
                        repParRBiax_NameB_temp = Node_parTrR{m,1}; % 1st
                        repParRBiax_NameB = repParRBiax_NameB_temp{n,1}(FormatY+1:end);          
                        % Ensure end value is part of nodenumber, if not, it is omitted
                        if repParRBiax_NameB == NodeNber{end,1}(FormatY+1:end)
                             n = n-1; % skip the specific iteration
                             continue;
                        end
                        repParRBiax{n, m} = str2num(strcat(repParRBiax_NameA,repParRBiax_NameB)); % 1st node, Element name
                        Elt_ParRBiax = [num2str(repParRBiax{n, m}), ', ', num2str(repParRBiax{n, m}), ', ' ,num2str(repParRBiax{n, m}+1) , ', ',  num2str(repParRBiax{n, m} +1 + dif_el), ', ', num2str(repParRBiax{n, m} + dif_el) ];
                        fprintf(fileID, '%s\r\n', Elt_ParRBiax);
                    end
                end
            end
            for m = 1:length(Node_hty)
                for n  = 1:length(Node_hty{m,1})
                    % Layer name
                    if comp_hty(m,1) == biax
                        htyBiax = cell2mat(Node_hty{m, 1}(n,1));
                        htyBiax_Name = htyBiax(FormatY+1:end);
                        if htyBiax_Name == NodeNber{end,1}(FormatY+1:end)
                             n = n-1; % skip the specific iteration
                             continue;
                        end
                        Elt_hty = [cell2mat(Node_hty{m, 1}(n,1)), ', ', cell2mat(Node_hty{m, 1}(n,1)), ', ' ,num2str(str2double(Node_hty{m,1}(n,1)) +1) , ', ',  num2str(str2double(Node_hty{m, 1}(n,1)) +1 + dif_el), ', ', num2str(str2double(Node_hty{m, 1}(n,1)) + dif_el) ];
                        fprintf(fileID, '%s\r\n', Elt_hty);
                    end
                end
            end            
            end
    end
end

%% BC: Points (to be changed later on for master/slave nodes
% BC nodes
Tlb = y_extTot{1, 1}(1)/1000; % Left edge bottom point, y coordinate
Trb = y_extTot{1, 1}(end)/1000; % Right edge bottom point, y coordinate
Tlt = y_extTot{1, LineNberExt}(1)/1000; % Left edge top point, y coordinate
Trt = y_extTot{1, LineNberExt}(end)/1000; % Right edge top point, y coordinate 
Lleft  = x_extTot{1, 1}(1)/1000; % Left edge point, x coordinate
Lright = x_extTot{1, 1}(end)/1000; % Right edge point, x coordinate

% (0,0)
PointNr1 = [num2str(1, NberFormaty), num2str(1, NberFormatx)];
Nr1X = num2str(Lleft, '%02.8f'); 
NrlY = num2str(Tlb, '%02.8f');

% (0,y_top)
PointNr2 = [num2str(LineNberExt, NberFormaty), num2str(1, NberFormatx)];
Nr2X = num2str(Lleft, '%02.8f'); 
Nr2Y = num2str(Tlt, '%02.8f');

% (x_end,y_top)
PointNr3 = [num2str(LineNberExt, NberFormaty), num2str(X_nodenber, NberFormatx)];
Nr3X = num2str(Lright, '%02.8f'); 
Nr3Y = num2str(Trt, '%02.8f');

% (x_end,0)
PointNr4 = [num2str(1, NberFormaty), num2str(X_nodenber, NberFormatx)];
Nr4X = num2str(Lright, '%02.8f'); 
Nr4Y = num2str(Trb, '%02.8f');

LineText = '**-------------------------------------';
fprintf(fileID, '%s\r\n', LineText); 
LineText =	'** PointNr1';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=PointNr1'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = PointNr1;
fprintf(fileID, '%s\r\n', LineText);
LineText = '**-------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** PointNr2';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=PointNr2'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = PointNr2;
fprintf(fileID, '%s\r\n', LineText);
LineText = '**-------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** PointNr3';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=PointNr3'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = PointNr3;
fprintf(fileID, '%s\r\n', LineText);
LineText =	'**------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText =	'** PointNr4';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=PointNr4'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = PointNr4;
fprintf(fileID, '%s\r\n', LineText); 

%% BC: Edges (Mid-point and Edges points)
% Left edge
SetId = x_mat == min(x_mat);
SetName = NodeNber(SetId, 1);

LineText = '**--------------------------------------'; 	
fprintf(fileID, '%s\r\n', LineText);
LineText = '**--------------------------------------'; 	 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Master node';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 

% Left Edge Points
if over_olap == 1
    LineNberEdge = LineNberExt -2 ; % subtract the added overlamina nodes
else
    LineNberEdge = LineNberExt;
end

LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=LeftEdgeNodeSet'; 
fprintf(fileID, '%s\r\n', LineText);
for i13 = 1 : LineNberEdge
    LineText = [];
    LineText = [LineText, SetName{i13}, ', '];
    fprintf(fileID, '%s\r\n', LineText);
end

% Right edge
SetId = x_mat == max(x_mat);
SetName = NodeNber(SetId, 1);

LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Master node';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 

% Right Edge Points
LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=RightEdgeNodeSet'; 
fprintf(fileID, '%s\r\n', LineText);
for ii13 = 1:LineNberEdge
    LineText = [];
    LineText = [LineText, SetName{ii13}, ', '];
    fprintf(fileID, '%s\r\n', LineText);
end

%% Sections and Orientations

if (mat_change == 1) && (flag==0)
    if over_lam == 1 && t_overmat == biax
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-8'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '** Section: Section-5-Eoverlam'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=BIAXMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        fprintf(fileID, '%s\r\n', LineText);        
    else if over_lam == 1 && t_overmat == ud
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-8'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '** Section: Section-5-Eoverlam'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=UDMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 
            fprintf(fileID, '%s\r\n', LineText);
        end       
    end
    
    LineText = '*Orientation, name=Ori-4'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-2-EPARENTBIAX'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Solid Section, elset=EPARENTBIAX, orientation=Ori-4, material=BIAXMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = ','; 
    fprintf(fileID, '%s\r\n', LineText);

    LineText = '*Orientation, name=Ori-2'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-5-EINTERFACE'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Cohesive Section, elset=EINTERFACE, orientation=Ori-2, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = strcat(num2str(t_Hres*(10^-3)),',');
    fprintf(fileID, '%s\r\n', LineText);

    LineText = '*Orientation, name=Ori-5'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-1-EPARENTUD'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Solid Section, elset=EPARENTUD, orientation=Ori-5, material=UDMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = ','; 

    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-3'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          0.,           1.,           0.,           1.,           0.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-6-ERESIN'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Cohesive Section, elset=ERESIN, orientation=Ori-3, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = strcat(num2str(t_Vres*(10^-3)),','); 

   
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-7'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-3-EREPAIRBIAX1'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Solid Section, elset=EREPAIRUD1, orientation=Ori-7, material=BIAXMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = ','; 

    
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-6'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '** Section: Section-4-EREPAIRUD1'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Solid Section, elset=EREPAIRBIAX1, orientation=Ori-6, material=UDMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = ','; 
    fprintf(fileID, '%s\r\n', LineText);
    
    LineText = '*Orientation, name=Ori-9'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '** Section: Section-4-Eresinrich'; 
    fprintf(fileID, '%s\r\n', LineText);
    %LineText = '*Cohesive Section, elset=Eresinrich, orientation=Ori-9, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
    LineText = '*Solid Section, elset=Eresinrich, orientation=Ori-9, material=RESINMATERIAL2'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = ','; 
    fprintf(fileID, '%s\r\n', LineText);
    
else if (mat_change == 1) && (flag==1)
    if over_lam == 1 && t_overmat == biax
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-8'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '** Section: Section-5-Eoverlam'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=BIAXMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        fprintf(fileID, '%s\r\n', LineText);        
    else if over_lam == 1 && t_overmat == ud
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-8'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '** Section: Section-5-Eoverlam'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=UDMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 
            fprintf(fileID, '%s\r\n', LineText);
        end       
    end
    
    LineText = '*Orientation, name=Ori-4'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-2-EPARENTBIAX'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Solid Section, elset=EPARENTBIAX, orientation=Ori-4, material=BIAXMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = ','; 
    fprintf(fileID, '%s\r\n', LineText);

    LineText = '*Orientation, name=Ori-2'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-5-EINTERFACE'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Cohesive Section, elset=EINTERFACE, orientation=Ori-2, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = strcat(num2str(t_Hres*(10^-3)),','); 
    fprintf(fileID, '%s\r\n', LineText);

    LineText = '*Orientation, name=Ori-5'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-1-EPARENTUD'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Solid Section, elset=EPARENTUD, orientation=Ori-5, material=UDMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = ','; 

    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-3'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          0.,           1.,           0.,           1.,           0.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-6-ERESIN'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Cohesive Section, elset=ERESIN, orientation=Ori-3, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = strcat(num2str(t_Vres*(10^-3)),','); 

    if layup(1,1) == biax
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-7'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-3-EREPAIRBIAX1'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Solid Section, elset=EREPAIRBIAX1, orientation=Ori-7, material=BIAXMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = ','; 
    end
    
    if layup(1,1) == ud
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-6'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '** Section: Section-4-EREPAIRUD1'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Solid Section, elset=EREPAIRUD1, orientation=Ori-6, material=UDMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = ','; 
    fprintf(fileID, '%s\r\n', LineText);
    end
    
    LineText = '*Orientation, name=Ori-9'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '** Section: Section-4-Eresinrich'; 
    fprintf(fileID, '%s\r\n', LineText);
    %LineText = '*Cohesive Section, elset=Eresinrich, orientation=Ori-9, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
    LineText = '*Solid Section, elset=Eresinrich, orientation=Ori-9, material=RESINMATERIAL2'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = ','; 
    fprintf(fileID, '%s\r\n', LineText);
    else if (flag==0)     
        if over_lam == 1 && t_overmat == biax
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-8'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '** Section: Section-5-Eoverlam'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=BIAXMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 
            fprintf(fileID, '%s\r\n', LineText);        
        else if over_lam == 1 && t_overmat == ud
                fprintf(fileID, '%s\r\n', LineText);
                LineText = '*Orientation, name=Ori-8'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '3, 0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '** Section: Section-5-Eoverlam'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=UDMATERIAL'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = ','; 
                fprintf(fileID, '%s\r\n', LineText);
            end       
        end
        
        LineText = '*Orientation, name=Ori-4'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** Section: Section-2-EPARENTBIAX'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=EPARENTBIAX, orientation=Ori-4, material=BIAXMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        fprintf(fileID, '%s\r\n', LineText);

        LineText = '*Orientation, name=Ori-2'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** Section: Section-5-EINTERFACE'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Cohesive Section, elset=EINTERFACE, orientation=Ori-2, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = strcat(num2str(t_Hres*(10^-3)),','); 
        fprintf(fileID, '%s\r\n', LineText);

        LineText = '*Orientation, name=Ori-5'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** Section: Section-1-EPARENTUD'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=EPARENTUD, orientation=Ori-5, material=UDMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 

        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-3'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          0.,           1.,           0.,           1.,           0.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** Section: Section-6-ERESIN'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Cohesive Section, elset=ERESIN, orientation=Ori-3, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = strcat(num2str(t_Vres*(10^-3)),','); 

        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-7'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** Section: Section-3-EREPAIRUD'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Solid Section, elset=EREPAIRUD, orientation=Ori-7, material=UDMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = ','; 

        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-6'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '** Section: Section-4-EREPAIRBIAX'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=EREPAIRBIAX, orientation=Ori-6, material=BIAXMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        fprintf(fileID, '%s\r\n', LineText);

    else if (flag==1) && layup(1,1) == ud
            if over_lam == 1 && t_overmat == biax
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-8'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '** Section: Section-5-Eoverlam'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=BIAXMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        fprintf(fileID, '%s\r\n', LineText);        
    else if over_lam == 1 && t_overmat == ud
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-8'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '** Section: Section-5-Eoverlam'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=UDMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 
            fprintf(fileID, '%s\r\n', LineText);
        end       
    end
            LineText = '*Orientation, name=Ori-2'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** Section: Section-5-EINTERFACE'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Cohesive Section, elset=EINTERFACE, orientation=Ori-2, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = strcat(num2str(t_Hres*(10^-3)),','); 
            fprintf(fileID, '%s\r\n', LineText);

            LineText = '*Orientation, name=Ori-5'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** Section: Section-1-EPARENTUD'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=EPARENTUD, orientation=Ori-5, material=UDMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 

            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-3'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          0.,           1.,           0.,           1.,           0.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** Section: Section-6-ERESIN'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Cohesive Section, elset=ERESIN, orientation=Ori-3, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = strcat(num2str(t_Vres*(10^-3)),','); 

            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-7'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** Section: Section-3-EREPAIRUD'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Solid Section, elset=EREPAIRUD, orientation=Ori-7, material=UDMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = ',';

        else if (flag==1) && layup(1,1) == biax
                if over_lam == 1 && t_overmat == biax
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-8'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '** Section: Section-5-Eoverlam'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=BIAXMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        fprintf(fileID, '%s\r\n', LineText);        
    else if over_lam == 1 && t_overmat == ud
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-8'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '** Section: Section-5-Eoverlam'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=Eoverlam, orientation=Ori-8, material=UDMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 
            fprintf(fileID, '%s\r\n', LineText);
        end       
    end
                LineText = '*Orientation, name=Ori-4'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '3, 0.'; 
                fprintf(fileID, '%s\r\n', LineText);
                LineText = '** Section: Section-2-EPARENTBIAX'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '*Solid Section, elset=EPARENTBIAX, orientation=Ori-4, material=BIAXMATERIAL'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = ','; 
                fprintf(fileID, '%s\r\n', LineText);

                LineText = '*Orientation, name=Ori-2'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '3, 0.'; 
                fprintf(fileID, '%s\r\n', LineText);
                LineText = '** Section: Section-5-EINTERFACE'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '*Cohesive Section, elset=EINTERFACE, orientation=Ori-2, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = strcat(num2str(t_Hres*(10^-3)),','); 
                fprintf(fileID, '%s\r\n', LineText);

                fprintf(fileID, '%s\r\n', LineText);
                LineText = '*Orientation, name=Ori-3'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '          0.,           1.,           0.,           1.,           0.,           0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '3, 0.'; 
                fprintf(fileID, '%s\r\n', LineText);
                LineText = '** Section: Section-6-ERESIN'; 
                fprintf(fileID, '%s\r\n', LineText);
                LineText = '*Cohesive Section, elset=ERESIN, orientation=Ori-3, controls=EC-1, material=RESIN, response=TRACTION SEPARATION, thickness=SPECIFIED'; 
                fprintf(fileID, '%s\r\n', LineText);
                LineText = strcat(num2str(t_Vres*(10^-3)),','); 

                fprintf(fileID, '%s\r\n', LineText);
                LineText = '*Orientation, name=Ori-6'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '3, 0.'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '** Section: Section-4-EREPAIRBIAX'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = '*Solid Section, elset=EREPAIRBIAX, orientation=Ori-6, material=BIAXMATERIAL'; 
                fprintf(fileID, '%s\r\n', LineText); 
                LineText = ','; 
                fprintf(fileID, '%s\r\n', LineText);
                LineText = ','; 
                fprintf(fileID, '%s\r\n', LineText);
            end
        end
    end  
    end
end


%% End Part
LineText = '*End Part'; 
fprintf(fileID, '%s\r\n', LineText);

%% Assembly and Instance 
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** ASSEMBLY'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Assembly, name=Assembly'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Instance, name=REPAIREDPANEL-1, part=REPAIREDPANEL';
fprintf(fileID, '%s\r\n', LineText);

%% End Instance
LineText = '*End Instance'; 
fprintf(fileID, '%s\r\n', LineText);

%% Node Set for BCs and BCs
LineText = '*Node'; 
fprintf(fileID, '%s\r\n', LineText);
LineTextRP_left = [num2str(1),', ', num2str(-0.001),', ', num2str((Tot_thick/2)/1000),', ',num2str(0.)  ];
fprintf(fileID, '%s\r\n', LineTextRP_left);
LineText = '*Node'; 
fprintf(fileID, '%s\r\n', LineText);
LineTextRP_right = [num2str(2),', ', num2str((Tot_length/1000)+0.001),', ', num2str((Tot_thick/2)/1000),', ',num2str(0.)  ];
fprintf(fileID, '%s\r\n', LineTextRP_right);
LineText = '*Nset, nset=RP_left'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '1,'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Nset, nset=RP_right'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '2,'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Constraint: ConstraintLeft'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Rigid Body, ref node=RP_left, tie nset=REPAIREDPANEL-1.LEFTEDGENODESET'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Constraint: ConstraintRight'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Rigid Body, ref node=RP_right, tie nset=REPAIREDPANEL-1.RIGHTEDGENODESET'; 
fprintf(fileID, '%s\r\n', LineText);

%% Crack Enrichment
if (flag==0) && (mat_change==0)
LineText = '*Enrichment, name=Crack-1, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EPARENTBIAX'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-2, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EPARENTUD'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-3, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRBIAX';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-4, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRUD'; 
fprintf(fileID, '%s\r\n', LineText);


else if(flag==0) && (mat_change==1)
LineText = '*Enrichment, name=Crack-1, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EPARENTBIAX'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-2, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EPARENTUD'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-3, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRBIAX1';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-4, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRUD1'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Enrichment, name=Crack-5, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.Eresinrich'; 
fprintf(fileID, '%s\r\n', LineText);
        
else if (flag==1) && layup(1,1) == ud
        LineText = '*Enrichment, name=Crack-2, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EPARENTUD'; 
        fprintf(fileID, '%s\r\n', LineText);
        if (mat_change == 0)
            LineText = '*Enrichment, name=Crack-4, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRUD';
            fprintf(fileID, '%s\r\n', LineText);
        end
        if (mat_change == 1)
            LineText = '*Enrichment, name=Crack-4, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRUD1';
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Enrichment, name=Crack-5, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.Eresinrich'; 
            fprintf(fileID, '%s\r\n', LineText);
        end
    else if (flag==1) && layup(1,1) == biax
            LineText = '*Enrichment, name=Crack-1, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EPARENTBIAX'; 
            fprintf(fileID, '%s\r\n', LineText);
            if (mat_change == 0)
                LineText = '*Enrichment, name=Crack-3, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRBIAX';
                fprintf(fileID, '%s\r\n', LineText);
            end
            if (mat_change == 1)
            LineText = '*Enrichment, name=Crack-4, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EREPAIRBIAX1';
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Enrichment, name=Crack-5, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.Eresinrich'; 
            fprintf(fileID, '%s\r\n', LineText);
            end
            end
        end
    end
    end


%% End of assembly
LineText = '*End Assembly'; 
fprintf(fileID, '%s\r\n', LineText);

%% Control
LineText = '** '; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** ELEMENT CONTROLS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** '; 
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat('*Section Controls, name=EC-1, ELEMENT DELETION=YES, KINEMATICS=AVERAGE STRAIN, VISCOSITY=',num2str(Viscosity)); 
%LineText = strcat('*Section Controls, name=EC-1, ELEMENT DELETION=YES, hourglass=STIFFNESS, VISCOSITY=',num2str(Viscosity)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = '1., 1., 1.'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** '; 
fprintf(fileID, '%s\r\n', LineText);

%% Materials
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Materials'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Material, name=UDMATERIAL'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Initiation, criterion=MAXPS';%Set CRITERION=MAXPS to specify a damage initiation criterion based on the maximum principal stress criterion for enriched elements. 
fprintf(fileID, '%s\r\n', LineText);
%LineText = ' 1e+09,'; %Maximum principal stress at damage initiation.
LineText = strcat(' ',num2str(Daminit_UD),','); %Maximum principal stress at damage initiation.
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat('*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=',num2str(Damev_UD)); 
fprintf(fileID, '%s\r\n', LineText);
%LineText = '100.,100.,100.'; 
LineText = strcat(num2str(Frac_UD1),',',num2str(Frac_UD2),',',num2str(Frac_UD3)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Elastic, type=ENGINEERING CONSTANTS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText =  strcat(num2str(E1UD),',',num2str(E2UD),',',num2str(E3UD),',',num2str(V12UD),',',num2str(V13UD),',',num2str(V23UD),',',num2str(G12UD),',',num2str(G13UD)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = num2str(G23UD);
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Material, name=BIAXMATERIAL'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Initiation, criterion=MAXPS'; %Set CRITERION=MAXPS to specify a damage initiation criterion based on the maximum principal stress criterion for enriched elements.
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat(' ',num2str(Daminit_BIAX),','); %Maximum principal stress at damage initiation.
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat('*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=',num2str(Damev_BIAX)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat(num2str(Frac_BIAX1),',',num2str(Frac_BIAX2),',',num2str(Frac_BIAX3));
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Elastic, type=ENGINEERING CONSTANTS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat(num2str(E1BIAX),',',num2str(E2BIAX),',',num2str(E3BIAX),',',num2str(V12BIAX),',',num2str(V13BIAX),',',num2str(V23BIAX),',',num2str(G12BIAX),',',num2str(G13BIAX)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = num2str(G23BIAX); 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Material, name=RESIN'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Initiation, criterion=QUADS'; %Set CRITERION=QUADS to specify a damage initiation based on the quadratic traction-interaction criterion for cohesive elements or enriched elements.
fprintf(fileID, '%s\r\n', LineText);
%LineText = ' 7.5e+06, 1.5e+07, 1.5e+07';%Maximum nominal stress in the normal-only mode. 
LineText = strcat(' ',num2str(Daminit_Resin1),',',' ',num2str(Daminit_Resin2),',',' ',num2str(Daminit_Resin3));%Maximum nominal stress in the normal-only mode. 
fprintf(fileID, '%s\r\n', LineText);
%Set TYPE=ENERGY to define the evolution of damage in terms of the energy required for failure (fracture energy) after the initiation of damage
%Set MIXED MODE BEHAVIOR=BK to specify the fracture energy as a function of the mode mix by means of the Benzeggagh-Kenane mixed mode fracture criterion.
%Set power this parameter equal to the exponent in the power law or the Benzeggagh-Kenane criterion that defines the variation of fracture energy with mode mix for cohesive elements.
LineText = strcat('*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=',num2str(Damev_Resin)); 
fprintf(fileID, '%s\r\n', LineText);
%Fracture energy.
LineText = strcat(num2str(Frac_Resin1),',',num2str(Frac_Resin2),',',num2str(Frac_Resin3));
fprintf(fileID, '%s\r\n', LineText);
%Damage stabilization
LineText = '*Damage Stabilization';
fprintf(fileID, '%s\r\n', LineText);
LineText = num2str(Damstab);
fprintf(fileID, '%s\r\n', LineText);
% *Elastic, type=TRACTION Use the following option to define uncoupled traction-separation behavior:
LineText = '*Elastic, type=TRACTION'; 
fprintf(fileID, '%s\r\n', LineText);
%Enn Ess and Ett for cohesive elements.
%LineText = ' 3e+14, 1.15e+14, 1.15e+14'; 
LineText = strcat(' ',num2str(Knn),',',' ',num2str(Kss),',',' ',num2str(Ktt));
fprintf(fileID, '%s\r\n', LineText);

LineText = '*Material, name=RESINMATERIAL2'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Initiation, criterion=MAXPS';%Set CRITERION=MAXPS to specify a damage initiation criterion based on the maximum principal stress criterion for enriched elements. 
fprintf(fileID, '%s\r\n', LineText);
%LineText = ' 1e+09,'; %Maximum principal stress at damage initiation.
LineText = strcat(' ',num2str(Daminit_Resin1),','); %Maximum principal stress at damage initiation.
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat('*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=',num2str(Damev_Resin)); 
fprintf(fileID, '%s\r\n', LineText);
%LineText = '100.,100.,100.'; 
LineText = strcat(num2str(Frac_Resin1),',',num2str(Frac_Resin2),',',num2str(Frac_Resin3)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Elastic, type=ENGINEERING CONSTANTS'; 
fprintf(fileID, '%s\r\n', LineText);
% I assume that the resin is an isotropic material and that the properties
% for cohesive resin and the material properties are the same.V12RESIN,E1RESIN
LineText =  strcat(num2str(E1RESIN),',',num2str(E1RESIN),',',num2str(E1RESIN),',',num2str(V12RESIN),',',num2str(V12RESIN),',',num2str(V12RESIN),',',num2str(E1RESIN/(2*(1+V12RESIN))),',',num2str(E1RESIN/(2*(1+V12RESIN)))); 
%LineText =  strcat(num2str(Knn),',',num2str(Knn),',',num2str(Knn),',',num2str(1/(Kss*2/Knn)-1),',',num2str(1/(Kss*2/Knn)-1),',',num2str(1/(Kss*2/Knn)-1),',',num2str(Kss),',',num2str(Ktt)); 
fprintf(fileID, '%s\r\n', LineText);
LineText = num2str(E1RESIN/(2*(1+V12RESIN)));
%LineText = num2str(Kss);
fprintf(fileID, '%s\r\n', LineText);

%% Step
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** STEP: Step-1'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Step, name=Step-1, nlgeom=YES, inc=10000'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Static'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat(num2str(first_Step),',',num2str(Tperiod),',',num2str(Min_incr),',',num2str(Max_incr));
fprintf(fileID, '%s\r\n', LineText);
LineText = '*CONTROLS, PARAMETERS=TIME INCREMENTATION';
fprintf(fileID, '%s\r\n', LineText);
LineText = strcat(num2str(I_o),',',num2str(I_R));
fprintf(fileID, '%s\r\n', LineText);

%% Boundary Condition
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** BOUNDARY CONDITIONS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Name: BCleft Type: Displacement/Rotation'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Boundary'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'RP_left, 1, 1'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'RP_left, 2, 2';
fprintf(fileID, '%s\r\n', LineText);
LineText = 'RP_left, 6, 6'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** Name: BCright Type: Displacement/Rotation'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Boundary'; 
fprintf(fileID, '%s\r\n', LineText);
LineText_Disp = ['RP_right' ,', ', num2str(1),', ', num2str(1),', ', num2str(x_disp/1000)];
fprintf(fileID, '%s\r\n', LineText_Disp); 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'RP_right, 2, 2'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'RP_right, 6, 6'; 
fprintf(fileID, '%s\r\n', LineText);

%% Field Output
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** OUTPUT REQUESTS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Restart, write, frequency=10'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** FIELD OUTPUT: F-Output-2'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Output, field'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Element Output, elset=REPAIREDPANEL-1.EINTERFACE, directions=YES'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'QUADSCRT, SDEG, STATUS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** FIELD OUTPUT: F-Output-3'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Element Output, elset=REPAIREDPANEL-1.ERESIN, directions=YES'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'QUADSCRT, SDEG, STATUS'; 
fprintf(fileID, '%s\r\n', LineText);

if mat_change == 1
%     fprintf(fileID, '%s\r\n', LineText);
%     LineText = '**'; 
%     fprintf(fileID, '%s\r\n', LineText);
%     LineText = '** FIELD OUTPUT: F-Output-3'; 
%     fprintf(fileID, '%s\r\n', LineText);
%     LineText = '**'; 
%     fprintf(fileID, '%s\r\n', LineText);
%     LineText = '*Element Output, elset=REPAIREDPANEL-1.Eresinrich, directions=YES'; 
%     fprintf(fileID, '%s\r\n', LineText);
%     LineText = 'QUADSCRT, SDEG, STATUS'; 
%     fprintf(fileID, '%s\r\n', LineText);
end

if (flag==0)
    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** FIELD OUTPUT: F-Output-4'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Element Output, elset=REPAIREDPANEL-1.EPARENTBIAX, directions=YES'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = 'STATUSXFEM, '; 
    fprintf(fileID, '%s\r\n', LineText);

    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** FIELD OUTPUT: F-Output-5'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Element Output, elset=REPAIREDPANEL-1.EPARENTUD, directions=YES'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = 'STATUSXFEM, '; 
    fprintf(fileID, '%s\r\n', LineText);
else if (flag==1) && layup(1,1) == ud
        LineText = '**'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** FIELD OUTPUT: F-Output-5'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '**'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Element Output, elset=REPAIREDPANEL-1.EPARENTUD, directions=YES'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = 'STATUSXFEM, '; 
        fprintf(fileID, '%s\r\n', LineText);
    else if (flag==1) && layup(1,1) == biax
            LineText = '**'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** FIELD OUTPUT: F-Output-4'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '**'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Element Output, elset=REPAIREDPANEL-1.EPARENTBIAX, directions=YES'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = 'STATUSXFEM, '; 
            fprintf(fileID, '%s\r\n', LineText);
        end
    end
end

LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** FIELD OUTPUT: F-Output-1'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Output, field, variable=PRESELECT'; 
fprintf(fileID, '%s\r\n', LineText);

%% History Output
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
% LineText = '** HISTORY OUTPUT: H-Output-1'; 
% fprintf(fileID, '%s\r\n', LineText);
% LineText = '**'; 
% fprintf(fileID, '%s\r\n', LineText);
% LineText = '*Output, history'; 
% fprintf(fileID, '%s\r\n', LineText);
% LineText = '*Node Output, nset=REPAIREDPANEL-1.RIGHTEDGECENTER'; 
% fprintf(fileID, '%s\r\n', LineText);
% LineText = 'RF1, U1'; 
% fprintf(fileID, '%s\r\n', LineText);
% LineText = '**'; 
% fprintf(fileID, '%s\r\n', LineText);
LineText = '** HISTORY OUTPUT: H-Output-1'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Output, history, variable=PRESELECT'; 
fprintf(fileID, '%s\r\n', LineText);

%% End Step
LineText = '*End Step'; 
fprintf(fileID, '%s\r\n', LineText);
end