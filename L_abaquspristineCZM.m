%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Load workspace from main
load('inputwithoutrepair.mat');

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

%% Sets: UD and BIAX
CompId_biax = zeros(TotElt, 1); 
CompId_ud = zeros(TotElt, 1); 
countElt = 1;
layup_el = layup(2:end); % bottom to top , bottom not taken since moved to y=0
Layup_El = (reshape([layup_el'; zeros(size(layup_el'))], 1, []))';

for i9 = 1 : LineNberExt - 1
    for i10 = 1 : X_nodenber - 1
        % Composite or not? 
        if mod(i9, 2) ~= 0 && Layup_El(i9) == biax
           CompId_biax(countElt, 1) = 2;  
        elseif mod(i9, 2) ~= 0 && Layup_El(i9) == ud
           CompId_ud(countElt, 1) = 3;
        end    
        countElt = countElt + 1;
    end
end    

if (flag==2)
    % BIAX
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Element, type=CPE4, Elset=EcompositeBiax';
    fprintf(fileID, '%s\r\n', LineText); 
    for ill = 1 : TotElt
        if CompId_biax(ill, 1) == 2
        LineText = EltName{ill, 1};
            for i12 = 1:size(EltNodes, 2) 
                LineText = [LineText,',',	EltNodes{ill, i12}];
            end
        fprintf(fileID, '%s\r\n', LineText);
        end
    end
    % UD
    LineText = '**--------------------------------------'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Element, type=CPE4, Elset=EcompositeUD';
    fprintf(fileID, '%s\r\n', LineText); 
    for ill = 1 : TotElt
        if CompId_ud(ill, 1) == 3
        LineText = EltName{ill, 1};
            for i12 = 1:size(EltNodes, 2) 
                LineText = [LineText,',',	EltNodes{ill, i12}];
            end
        fprintf(fileID, '%s\r\n', LineText);
        end
    end
else if (flag==3) && layup(1,1) == ud 
        % UD
        LineText = '**--------------------------------------'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Element, type=CPE4, Elset=EcompositeUD';
        fprintf(fileID, '%s\r\n', LineText); 
        for ill = 1 : TotElt
            if CompId_ud(ill, 1) == 3
            LineText = EltName{ill, 1};
                for i12 = 1:size(EltNodes, 2) 
                    LineText = [LineText,',',	EltNodes{ill, i12}];
                end
            fprintf(fileID, '%s\r\n', LineText);
            end
        end

    else if (flag==3) && layup(1,1) == biax
            % BIAX
            LineText = '**--------------------------------------'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Element, type=CPE4, Elset=EcompositeBiax';
            fprintf(fileID, '%s\r\n', LineText); 
            for ill = 1 : TotElt
                if CompId_biax(ill, 1) == 2
                LineText = EltName{ill, 1};
                    for i12 = 1:size(EltNodes, 2) 
                        LineText = [LineText,',',	EltNodes{ill, i12}];
                    end
                fprintf(fileID, '%s\r\n', LineText);
                end
            end
        end
    end
end

%% Sets: Resin Interface
LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Element,type=COH2D4, Elset=Einterface'; % Cohesive elements
fprintf(fileID, '%s\r\n', LineText); 
for ii11 = 1 : TotElt
    if CompId(ii11, 1) == 0
        LineText = EltName{ii11, 1};
            for ii12 = 1 : size(EltNodes, 2)
                LineText = [LineText,', ',	EltNodes{ii11, ii12}];
            end
        fprintf(fileID, '%s\r\n', LineText);
    end
end

%% BC: Points (to be changed later on for master/slave nodes
% BC nodes
Tlb = y_extTot{1, 1}(1); % Left edge bottom point, y coordinate
Trb = y_extTot{1, 1}(end); % Right edge bottom point, y coordinate
Tlt = y_extTot{1, LineNberExt}(1); % Left edge top point, y coordinate
Trt = y_extTot{1, LineNberExt}(end); % Right edge top point, y coordinate 
Lleft  = x_extTot{1, 1}(1); % Left edge point, x coordinate
Lright = x_extTot{1, 1}(end); % Right edge point, x coordinate

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
LineText = '**--------------------------------------'; 
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Node';
fprintf(fileID, '%s\r\n', LineText); 
LineText = '*Nset, nset=LeftEdgeNodeSet'; 
fprintf(fileID, '%s\r\n', LineText);
for i13 = 1 : LineNberExt
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
for ii13 = 1:LineNberExt
    LineText = [];
    LineText = [LineText, SetName{ii13}, ', '];
    fprintf(fileID, '%s\r\n', LineText);
end

%% Sections and Orientations
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
LineText = '0.08,'; 
fprintf(fileID, '%s\r\n', LineText); 

if (flag==2)
    LineText = '*Orientation, name=Ori-4'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-2-EcompositeBiax'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Solid Section, elset=EcompositeBiax, orientation=Ori-4, material=BIAXMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = ','; 
    
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-5'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** Section: Section-1-EcompositeUD'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '*Solid Section, elset=EcompositeUD, orientation=Ori-5, material=UDMATERIAL'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = ',';

    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Orientation, name=Ori-3'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
    fprintf(fileID, '%s\r\n', LineText); 
    LineText = '3, 0.'; 
    fprintf(fileID, '%s\r\n', LineText);
else if (flag==3) && layup(1,1) == ud
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-5'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** Section: Section-1-EcompositeUD'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '*Solid Section, elset=EcompositeUD, orientation=Ori-5, material=UDMATERIAL'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = ','; 
        
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Orientation, name=Ori-3'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
        fprintf(fileID, '%s\r\n', LineText); 
        LineText = '3, 0.'; 
        fprintf(fileID, '%s\r\n', LineText);
        
    else if (flag==3) && layup(1,1) == biax
            LineText = '*Orientation, name=Ori-4'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** Section: Section-2-EcompositeBiax'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '*Solid Section, elset=EcompositeBiax, orientation=Ori-4, material=BIAXMATERIAL'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = ','; 
            
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Orientation, name=Ori-3'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '          1.,           0.,           0.,           0.,           1.,           0.'; 
            fprintf(fileID, '%s\r\n', LineText); 
            LineText = '3, 0.'; 
            fprintf(fileID, '%s\r\n', LineText);
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

%% Constraint - Rigid Body
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
if (flag==2)
    LineText = '*Enrichment, name=Crack-1, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EcompositeUD'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Enrichment, name=Crack-2, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EcompositeBIAX'; 
    fprintf(fileID, '%s\r\n', LineText);
else if (flag==3) && layup(1,1) == ud
        LineText = '*Enrichment, name=Crack-1, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EcompositeUD'; 
        fprintf(fileID, '%s\r\n', LineText);
    else if (flag==3) && layup(1,1) == biax
            LineText = '*Enrichment, name=Crack-2, type=PROPAGATION CRACK, elset=REPAIREDPANEL-1.EcompositeBIAX'; 
            fprintf(fileID, '%s\r\n', LineText);
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
LineText = '*Section Controls, name=EC-1, ELEMENT DELETION=YES, KINEMATICS=AVERAGE STRAIN, VISCOSITY=1e-05'; 
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
LineText = '*Damage Initiation, criterion=MAXPS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 1e+09,'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=2.'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '100.,100.,100.'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Elastic, type=ENGINEERING CONSTANTS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 4.163e+10, 1.493e+10, 1.343e+10,     0.241,    0.2675,    0.3301, 5.047e+09, 5.047e+09'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 5.047e+09,';
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Material, name=BIAXMATERIAL'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Initiation, criterion=MAXPS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 1e+09,'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=2.'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '100.,100.,100.'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Elastic, type=ENGINEERING CONSTANTS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 1.392e+10, 1.392e+10, 1.392e+10,     0.53,    0.19,    0.02, 1.15e+10, 1.15e+10'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 4.539e+9,'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Material, name=RESIN'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Initiation, criterion=QUADS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 7.5e+06, 1.5e+07, 1.5e+07'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Damage Evolution, type=ENERGY, mixed mode behavior=BK, power=3.8'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '200.,1000.,1000.'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Elastic, type=TRACTION'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = ' 3e+14, 1.15e+14, 1.15e+14'; 
fprintf(fileID, '%s\r\n', LineText);

%% Step
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** STEP: Step-1'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Step, name=Step-1, nlgeom=YES, inc=100'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Static'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '1.0E-1,1.0,1.0E-15,1.0E-1'; 
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
LineText_thetaL = ['RP_left' ,', ', num2str(6),', ', num2str(6),', ', num2str(-theta_rot)];
fprintf(fileID, '%s\r\n', LineText_thetaL); 

LineText = '** Name: BCright Type: Displacement/Rotation'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Boundary'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = 'RP_right, 2, 2'; 
fprintf(fileID, '%s\r\n', LineText);
LineText_thetaR = ['RP_right' ,', ', num2str(6),', ', num2str(6),', ', num2str(theta_rot)];
fprintf(fileID, '%s\r\n', LineText_thetaR);

%% Field Output
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '** OUTPUT REQUESTS'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '**'; 
fprintf(fileID, '%s\r\n', LineText);
LineText = '*Restart, write, frequency=1'; 
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

if (flag==2)
    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** FIELD OUTPUT: F-Output-4'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Element Output, elset=REPAIREDPANEL-1.EcompositeBIAX, directions=YES'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = 'STATUSXFEM, '; 
    fprintf(fileID, '%s\r\n', LineText);

    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '** FIELD OUTPUT: F-Output-4'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '**'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = '*Element Output, elset=REPAIREDPANEL-1.EcompositeUD, directions=YES'; 
    fprintf(fileID, '%s\r\n', LineText);
    LineText = 'STATUSXFEM, '; 
    fprintf(fileID, '%s\r\n', LineText);
    
else if (flag==3) && layup(1,1) == ud
        LineText = '**'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '** FIELD OUTPUT: F-Output-4'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '**'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = '*Element Output, elset=REPAIREDPANEL-1.EcompositeUD, directions=YES'; 
        fprintf(fileID, '%s\r\n', LineText);
        LineText = 'STATUSXFEM, '; 
        fprintf(fileID, '%s\r\n', LineText);   
    else if (flag==3) && layup(1,1) == biax
            LineText = '**'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '** FIELD OUTPUT: F-Output-4'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '**'; 
            fprintf(fileID, '%s\r\n', LineText);
            LineText = '*Element Output, elset=REPAIREDPANEL-1.EcompositeBIAX, directions=YES'; 
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

%% Stiffness Calculation

% Create table for X_tot (r), Y_Tot (l0 laminae layers 9 layers of resin - bottom to top)
Ones = ones(length(Y_Tot(2,:)'),1);
columnLayer = zeros(length(Ones), LineNberExt-1);
cnt_lin = linspace(1,19,19);

for i = 1:length(layup)-1
    id = cnt_lin(1:2:end);
    columnLayer(:,id(i)) = (layup(i+1,:)-t_Hres).*Ones/1000;
end
for k = 1: length(layup)-2
    idx = cnt_lin(2:2:end);
    columnLayer(:,idx(k)) = t_Hres.*Ones/1000;
end

T = table(X_Tot(1,:)'/1000, columnLayer);
writetable(T,'ModelInput.xlsx') 

% Run stiffness calculation for no repair panel
run('StiffnessMain.m');
        

