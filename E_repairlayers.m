%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: Assemble repair layers

function [x_rep, y_rep, Y_top] = E_repairlayers(yrep_left_aseb, xrep_left_aseb, yrep_right_aseb, xrep_right_aseb, x_hty, y_right_par, y_left_par, L_remrow)
        for i9 = 1 : length(L_remrow)
            % Take end edge values of the repairs
            y_left_in(i9,:)   = yrep_left_aseb{i9,1}(1,1);
            x_left_in(i9,:)   = xrep_left_aseb{i9,1}(1,1);
            y_right_end(i9,:) = yrep_right_aseb{i9,1}(1,end);
            x_right_end(i9,:) = xrep_right_aseb{i9,1}(1,end);

            % Find index values on x_hty and take values from hty layers to build repair layers x_rep
            idx_rep_left(i9,:)  = find(x_left_in(i9,:)   == x_hty(1,:));
            idx_rep_right(i9,:) = find(x_right_end(i9,:) == x_hty(1,:));

            % Construct x_rep by taking coresponding elements at x_hty (bottom to top)
            x_rep{i9,1} = x_hty(1,idx_rep_left(i9,:)+1:idx_rep_right(i9,:)-1);
            
            % y_rep is created using linspace, since the corresponding y_hty value can be taken (different!!)    
            y_rep{i9,1} = linspace(y_left_in(i9,:), y_right_end(i9,:), length(x_rep{i9,1}));   
            
            % Concatenate parent and rep layers
            Y_top(i9,:) = [y_left_par{i9,1}(:,:)  y_rep{i9,1}(:,:)  y_right_par{i9,1}(:,:) ];    
        end
end