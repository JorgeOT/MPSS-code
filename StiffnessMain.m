%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)
Stiffnessmain_function =
%% Import Data
load inputwithoutrepair.mat
filename1 = 'ModelInput.xlsx'; % Top to Bottom layer stacking
filename2 = 'MaterialInput.xlsx';  
Layup = xlsread(filename1); % Top to Bottom layer stacking
material = xlsread(filename2);

r = Layup(:,1);         % [m]

%% Material Properties for the Biax and Uniax
% Uniax
E1_uniax  = material(1,1);           % Pa
E2_uniax  = material(2,1);           
E3_uniax  = material(3,1);          
G12_uniax = material(4,1);         
G13_uniax = material(5,1);          
G23_uniax = material(6,1);       
v12_uniax = material(7,1);           % [-]
v21_uniax = material(8,1);           % [-]
v13_uniax = material(9,1);           % [-]

% Biax
E1_biax   = material(1,2);           % [N/m^2]  Pa
E2_biax   = material(2,2);          
E3_biax   = material(3,2);        
G12_biax  = material(4,2);       
G13_biax  = material(5,2);        
G23_biax  = material(6,2);         
v12_biax  = material(7,2);    
v21_biax  = material(8,2);           
v13_biax  = material(9,2);                  

% Resin
E1_resin  = material(1,3);            % Pa
E2_resin  = material(2,3); 
G12_resin = material(4,3);
v12_resin = material(7,3);
v21_resin = material(8,3);

%% Transformation matrix---------------------------------------------------
theta = 0;                            % In Degrees

c = cos(theta);
s = sin(theta);
 
T = [c^2 s^2 -2*s*c; 
     s^2 c^2  2*s*c;
     s*c -s*c c^2-s^2];
 
Transpose_T = [c^2 s^2 s*c; 
               s^2 c^2  -s*c;
               -2*s*c 2*s*c c^2-s^2];
           
% Local stiffnes Matrix for the Uniax and Biax-----------------------------
Ql_uniax = (1/(1-v12_uniax*v21_uniax))*[E1_uniax  v21_uniax*E1_uniax  0;
                                        v12_uniax*E2_uniax  E2_uniax  0;
                                        0 0 G12_uniax*(1-(v12_uniax*v21_uniax))];

Ql_biax  = (1/(1-v12_biax*v21_biax))*[E1_biax  (v21_biax*E1_biax)  0;
                                      (v12_biax*E2_biax)  E2_biax  0;
                                      0 0 (G12_biax*(1-(v12_biax*v21_biax)))];
                                  
Ql_resin = (1/(1-v12_resin*v21_resin))*[E1_resin  (v21_resin*E1_resin)  0;
                                      (v12_resin*E2_resin)  E2_resin  0;
                                      0 0 (G12_resin*(1-(v12_resin*v21_resin)))];

% Global stifness of uniax and biax----------------------------------------
Qg_uniax = T* Ql_uniax *Transpose_T;
Qg_biax  = T* Ql_biax  *Transpose_T;
Qg_resin = T* Ql_resin *Transpose_T;

% Thicknesses of each layer at each point of the blade---------------------

for i = 1:length(Layup(1,:))-1 % Top to Bottom layer stacking
    h(:,i) = Layup(:,i+1); % [m]   
end

% Total thickness
t = sum(h,2); % [m]   

% Number of plies and resin layers
Totply = length(Y_res_Tot(:,1)) - 1;
Totres = Totply - 1;
Totlay = Totply + Totres;

if ~mod(Totply,2) % Even Number of plies
    for i = 1: Totply
        z_top(:,i) = h(:,i);
    end
    z_top(:,end) = z_top(:,end)/2;
    Z_top = sum(z_top,2).*ones(size(z_top));
    
    for i = 1:Totply-1
       Z_top(:,i+1) = Z_top(:,i) - z_top(:,i);
       Z_bot = -flip(Z_top,2);
    end
    Z = [Z_top Z_bot]; % [m]
end

[q,w] = size(Z);                     

for i = 1:q
    % For each layer distance with respect to neutral mid-axis
    for j = 1 : w-1
        diffA(i,j) = (Z(i,j)-Z(i,j+1));
        diffB(i,j) = (Z(i,j)^2-Z(i,j+1)^2);
        diffD(i,j) = ((Z(i,j)^3)-(Z(i,j+1)^3));
    end
    
    for m = 1:length(diffA(1,:))  
        if h(:,m) == t_biax/1000
            A_temp(:,:,m) = Qg_biax.*diffA(i,m);
            B_temp(:,:,m) = Qg_biax.*diffB(i,m);
            D_temp(:,:,m) = Qg_biax.*diffD(i,m);
        else if h(:,m) == t_ud/1000
                A_temp(:,:,m) = Qg_uniax.*diffA(i,m);
                B_temp(:,:,m) = Qg_uniax.*diffB(i,m);
                D_temp(:,:,m) = Qg_uniax.*diffD(i,m);
            else if h(:,m) == t_Hres/1000
                    A_temp(:,:,m) = Qg_resin.*diffA(i,m);
                    B_temp(:,:,m) = Qg_resin.*diffB(i,m);
                    D_temp(:,:,m) = Qg_resin.*diffD(i,m);
                end
            end
        end
    end
    
    % Composing the ABD Matrix
    A = sum(A_temp,3);
    B = (1/2)*sum(B_temp,3);
    D = (1/3)*sum(D_temp,3);   
       
    [ABD] = [A B;
             B D];
         
    invD = inv(D);
    
    % Calculate Ex since the layup is symetric.
    Ex   = (1./t).*(A(1,1)+(A(1,2)*((A(2,3)*A(1,3)-A(1,2)*A(3,3))/(A(2,2)*A(3,3)-A(2,3)^2))  + A(1,3)*((-A(1,3)/A(3,3)) + ((A(2,3)*A(1,2)*A(3,3))-A(2,3)^2*A(1,3))/((A(2,2)*A(3,3)^2-A(2,3)^2*A(3,3))))));
    Ey   = (1./t).*(A(2,2)+(A(1,2)*((A(2,3)*A(1,3)-A(1,2)*A(3,3))/(A(1,1)*A(3,3)-A(1,3)^2))  + A(2,3)*((-A(2,3)/A(3,3)) + ((A(1,3)*A(1,2)*A(3,3))-A(1,3)^2*A(2,3))/((A(1,1)*A(3,3)^2-A(1,3)^2*A(3,3))))));
    Gxy  = (1./t).*(A(3,3) - (A(2,3)^2/A(2,2)) + (((2*A(1,3)*A(1,2)*A(2,2)*A(2,3)) - ((A(1,2)^2)*(A(2,3)^2)) - ((A(1,3)^2) * (A(2,2)^2)))/ ((A(1,1)*A(2,2)^2) - A(1,2)^2*A(2,2))));
    Nuxy = (A(1,2) - ((A(1,3)*A(2,3))/A(3,3)))/(A(2,2) - (A(2,3)^2/A(3,3)));
    
    % Moment of inertia
    Ix = (r(end)*t(end)^3)/12;
    K  = Ex.* Ix;
    
    %% Flexural Rigidity
    
    if (flag == 2)
        % Number of biax and ud composite
        Totbiax = length(find(h(1,:) == (t_biax/1000)));
        Totud   = length(find(h(1,:) == (t_ud/1000)));
        
        % Flexural rigidity each layer component
        D_biax  = Totbiax * E1_biax  * (1/12 * r(end)* (t_biax/1000)^3); % [N/m^2 * m^4 = N.m^2]
        D_uniax = Totud   * E1_uniax * (1/12 * r(end)* (t_ud/1000)^3);
        D_resin = Totres  * E1_resin * (1/12 * r(end)* (t_Hres/1000)^3);
        D_add = D_biax + D_uniax + D_resin; % [N.m^2]
        
    else if (flag == 3) && layup(1,1) == ud
            % Number of biax and ud composite
            Totud   = length(find(h(1,:) == (t_ud/1000)));

            % Flexural rigidity each layer component
            D_uniax = Totud   * E1_uniax * (1/12 * r(end)* (t_ud/1000)^3); % [N/m^2 * m^4 = N.m^2]
            D_resin = Totres  * E1_resin * (1/12 * r(end)* (t_Hres/1000)^3);
            D_add = D_uniax + D_resin; % [N.m^2]
            
        else if (flag == 3) && layup(1,1) == biax
                % Number of biax and ud composite
                Totbiax = length(find(h(1,:) == (t_biax/1000))); 

                % Flexural rigidity each layer component
                D_biax  = Totbiax * E1_biax  * (1/12 * r(end)* (t_biax/1000)^3); % [N/m^2 * m^4 = N.m^2]
                D_resin = Totres  * E1_resin * (1/12 * r(end)* (t_Hres/1000)^3);
                D_add = D_biax + D_resin; % [N.m^2]
                
            end
        end
    end
    
    % Take the half of composite
    H = flip(h(:,1:(length(Y_Tot(:,1))/2)-1),2);
    distance1 = flip(h(:,1:(length(Y_Tot(:,1))/2)-1),2);                    % For storing inside the loop
    
    for j = 1: (length(Y_Tot(:,1))/2)-1
       mid = h(:,(length(Y_Tot(:,1))/2));                                   % Initial distance from symmetriy axis
       distance0(:,j) = H(:,j) + mid;                                       % Initial distance plus each layer of stacking (full thickness since we are taking half)
       distance1(:,j) = distance1(:,j)*2;                                   % Each layer is multiplied by 2 in order to take the other half of composite
    end

    for j = 1: (length(Y_Tot(:,1))/2)-2
        distance1(:,j+1) = distance1(:,j+1) + distance1(:,j);               % Accumulate the sum by each layer and store
        d(:,j) = distance0(:,j+1) + distance1(:,j);                         % Add distances
    end
    
    % Distances from symmetric neutral axis
    distance = [distance0(:,1)  d(:,:)]; % [m]                              % Construct the final distances
    
    % Calculate the flexural rigidity component of layers with respect to symmetric neautral axis
    for j = 1:(length(Y_Tot(:,1))/2)-1
        if H(:,j) == t_biax/1000
            Dflex_temp(:,j) = 1/2 * E1_biax * r(end)* (t_biax/1000) * distance(1,j)^2; % [N/m^2 * m * m^2 = N.m^2]
        else if H(:,j) == t_ud/1000
                Dflex_temp(:,j) = 1/2 * E1_uniax * r(end)* (t_ud/1000) * distance(1,j)^2;
            else if H(:,j) == t_Hres/1000
                    Dflex_temp(:,j) = 1/2 * E1_resin * r(end)* (t_Hres/1000) * distance(1,j)^2;
                end
            end
        end
    end
    
    % Flexural rigidity
    Dflex = sum(Dflex_temp) + D_add; % [N.m^2]
    
end

% Print values
formatSpec = 'CASE1A: Non-repaired panel';
fprintf(formatSpec);

Ex_val = Ex(1,1);
formatSpec = 'Analytical Young Modulus Ex under tension is %4.2f Pa \n';
fprintf(formatSpec,Ex_val);

Ey_val = Ey(1,1);
formatSpec = 'Young Modulus Ey is %4.2f Pa \n';
fprintf(formatSpec,Ey_val);

Gxy_val = Gxy(1,1);
formatSpec = 'Shear Modulus Gxy is %4.2f Pa \n';
fprintf(formatSpec,Gxy_val);

Nuxy_val = Nuxy(1,1);
formatSpec = 'Poisson ratio Nuxy is %4.2f [-] \n';
fprintf(formatSpec,Nuxy_val);

D_tot_val = Dflex;
formatSpec = 'Analytical Flexural Rigidity is %4.2f [Nm] \n';
fprintf(formatSpec,D_tot_val);


   
