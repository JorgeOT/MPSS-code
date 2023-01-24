%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: Take the last data point (y_start) from the bottom layer of the layers with defect 
% and use it as value reference to extend the healthy layers (dashed lines) by knowing the layup of the healthy 
% layers and the thicknesses of the layer material (contained in temps_lay that is subtracted to the matrix 
% created based on y_start).

function [x_aseb, y_aseb, x_hty, y_hty] = C_healthylayers(x_store, y_store, layup_nodefect, layup_defect)

        % Take last bottom defected layer and create a matrix
        y_start     = y_store(1,1);
        yline       = y_start * ones(length(layup_nodefect), length(y_store));

        % Prepare a matrix which contains information of the layup and thicknesses of lamina
        temp_nd_lay = layup_nodefect.* ones(length(layup_nodefect), length(yline)); 
        temp_line   = layup_defect(end)*ones(1, length(yline));
        temps_lay   = [temp_line; temp_nd_lay];

        % Subtract
        for i3 = 1 : length(layup_nodefect)
            temps_lay(i3 + 1,:) = temps_lay(i3 + 1,:) +  temps_lay(i3 ,:);
            YLines(i3,:)        = yline(i3,:) - temps_lay(i3,:);      
        end

        % Store healthy layers after extension and subtraction of layup information
        YLine = flip(YLines);
        XLine = repmat(x_store(1, :),length(layup_nodefect),1);

% Later on, all layers are moved to bottom (red arrow), to ensure bottom layer is at y = 0. 
% Now, we have y_aseb, x_aseb (assembled and moved down) and x_hty and y_hty (healthy layers).        
        
        % Assemble defect layers (x and y store) and healthy layers to include the bottom and top layer
        x_aseb = [ XLine; x_store]; 
        y_asebs = [ YLine; y_store];
        y_aseb = y_asebs - y_asebs(1,:);                                    % NOTE: Move y layers such that bottom layers y = 0

        % Assemble healthy laminates where bottom layers are moved to y=0
        x_hty = XLine; 
        y_hty = YLine- y_asebs(1,:);
end