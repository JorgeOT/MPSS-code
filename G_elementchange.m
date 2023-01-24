%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: This function create the coarse to fine elements, finding
% initially the mid x-values in each section (except the first and last
% sections since we only want coarse to fine mesh there instead of fine-coarse-fine for the remaining sections.)
% By knowing initial, mid and end value of each section, a logarithmic
% scale x-axis is created per half-section.

function [X_res_Tot, Y_res_Tot] = G_elementchange(width_vec, Tot_layup, xLsec, xRsec, xRsecEnd, Xad_left, Xad_right, Y_res_tot)
        % Concatenate all sections without resin (isolated)
        Xres_Tot_el = [xLsec xRsec xRsecEnd];
        
        % Take corresponding initial value wrt to resins
        x_elem_in = [ Xres_Tot_el{1,1}(:,1) Xad_left Xad_right];
        
            Width_vec = width_vec + 1; % +1 since codes works by taking for example: 5 x-points making 4 elements
            for i10 = 2: length(Width_vec)
                % For every cell take end value 
                x_elem_end(1,1)    = Xres_Tot_el{1,1}(1,end);
                x_elem_end(i10,1)  = Xres_Tot_el{1,i10}(1,end);
            end
            
            % Find middle value of all sections (except first and last section)
            for i11 = 1: length(Width_vec)
                x_elem{i11,1} = repmat(linspace(x_elem_in(1,i11),x_elem_end(i11,1), Width_vec(1,i11)),length(Tot_layup),1);
                
                % Mid value
                x_elem_mid(i11,1) = ((x_elem{i11,1}(1,end) - x_elem{i11,1}(1,1))/2) + x_elem{i11,1}(1,1);      
            end
            
            % Construct in, mid and end vlaues
            x_elem_edge = [x_elem_in(1,:)' x_elem_mid x_elem_end ];
            
            % For mid-section, take mid values at each sec
            x_elem_midsec = x_elem_edge(2:end-1,:); % not considering first and end sections
            w_vec = Width_vec(1,2:end-1)'; % not considering first and end sections
            
            X_Midsec = [];
            for i = 1:length(x_elem_midsec)
                % The first column (x_elem_in) is subtracted by 1, so that
                % when we make exp(linspace(log))) space, we start from 1
                % (midsec(1,1)) until value of 2nd column(x_elem_mid) subtracted by
                % diffmidsec. This is done since values do not start from
                % 1, making function not give good element meshing. Also,
                % to ensure that mid value is not doubled.
                
                x_elem_In(i,1) = x_elem_midsec(i,1) - 1;
                midsec(i,1) = 1; 
                midsec(i,2) = x_elem_midsec(i,2) - x_elem_In(i,1); 
                
                % Function that makes fine to coarse meshes
                incval_midsec{i,1} = exp(linspace(log(midsec(i,1)),log(midsec(i,2)),w_vec(i,1)));
                
                % Right part of mid-section: Take the vaues from function and subtract them in 3rd
                % column (x_elem_end) and add 1 (since a subraction of 1 was
                % done initially)
                R_midsec{i,1} = x_elem_midsec(i,3) - flip(incval_midsec{i,1}) + 1;
                
                % For the left part of the mid-section
                L_midsec{i,1} = incval_midsec{i,1} + x_elem_In(i,1);
                
                % concatenate left and right half sections
                Midsec{i,1} = repmat([L_midsec{i,1} R_midsec{i,1}(2:end)],length(Tot_layup),1); % R_midsec starts at element 2 since we do not want two same values in one array
                
                % Put into cell and concatetate them all
                X_Midsec = [X_Midsec Midsec{i,1}];
            end 
            % For the first and last section:
            
            % Take end width
            W_vec = [Width_vec(1,1) Width_vec(1,end)]';
            
            % Vary elements at end secs
            x_elem_ends = [x_elem_edge(1,:); x_elem_edge(end,:)];
            incval_endsec = exp(linspace(log(x_elem_ends(1,1)),log(x_elem_ends(1,end)),W_vec(1,1)));
            x_elem_end1 = repmat(x_elem_ends(1,end) - flip(incval_endsec) + 1,length(Tot_layup),1);
            x_elem_end2 = repmat(incval_endsec + x_elem_ends(2,1),length(Tot_layup),1)-1;
            
            % Concatenate for end secs
            X_res_Tot = [x_elem_end1 X_Midsec x_elem_end2];
            X_res_Tot = X_res_Tot - X_res_Tot(:,1);
            Y_res_Tot = Y_res_tot(:, 1:length(X_res_Tot));
end
        
