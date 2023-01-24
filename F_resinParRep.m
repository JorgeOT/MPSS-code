%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: Creates resin in between parent and laminate composites, by
% adding user defined thickness of the resin. Note that when the a resin is
% added, then consecutive values will be moved towards positive
% x-direction. Different sections 'sec' are created. First section
% represent the parent laminate at left panel. Section represent
% the section after the first resin addition, and so on (see Figure). This
% is all done until last section, achieving new x and y points with the
% addition of resin in between parent and repair.

function [xLsec, xRsec, xRsecEnd, Xad_left, Xad_right] = F_resinParRep(L_remrow, x_left_par, y_left_par, x_right_par, y_right_par, t_Vres, Tot_layup, x_rep, y_rep)

        xLsec = cell(1,length(L_remrow));
        xRsec = cell(1,length(L_remrow));
        
        % Add resin thickness on edge points of parent
        for ii9 = 1:length(L_remrow)
            
        % LEFT (Top to bottom)
            % Take end-value at each parent layer
            xval_res_left(ii9,1)  = x_left_par{end-ii9+1,1}(1,end);
            yval_res_left(ii9,1)  = y_left_par{end-ii9+1,1}(1,end);
  
            % Create resin elements for left panel by adding t_Vres thickness
            xres_left(ii9,1)  = xval_res_left(ii9,1)  + (t_Vres*ii9); % multiplied since all other values moved by t_Vres towards right
            yres_left(ii9,1)  = yval_res_left(ii9,1);

            % Build the x-axis which contains all resin in all layup (using Tot_layup)
            xad_left{1,ii9}   = repmat(xres_left(ii9,1),length(Tot_layup),1);
            yad_left{1,ii9}   = repmat(yres_left(ii9,1),length(Tot_layup),1);
            
            % Length of previous layer and +1 to take next element value in next laminate
            length_Lsec(ii9,1) = length(x_left_par{end-ii9+1,1}(:,:)) + 1;   % this is needed for later   

        % RIGHT (Bottom to top)
      
            % Take end-value at each right parent layer
            xval_res_right(ii9,1) = x_rep{ii9,1}(1,end);
            yval_res_right(ii9,1) = y_rep{ii9,1}(1,end);
            
            % Create resin in right panel by adding t_Vres thickness
            xres_right(ii9,1) = xval_res_right(ii9,1) + (t_Vres*ii9); % multiplied since all other values moved by t_Vres towards right
            yres_right(ii9,1) = yval_res_right(ii9,1);
            
            % Take initital value of right parent in order to locate it later
            idx_Rsec_temp(ii9,1) = x_right_par{ii9,:}(1,1);  % this is needed for later   
        end
        
        % Take initial values of the adehsives in the right panel then add resin due to moving of x values to right direction
        for ii9 = 1:length(L_remrow)
            xRrep_par(ii9,1) = xres_right(ii9,1) + (t_Vres*(length(L_remrow)-1));
        end
        
        for iii9 = 1:length(L_remrow)-1
            % From length information, we know where we take from which value and to which value in the next sections
            x_left_parR{iii9,1} = x_left_par{end-iii9,1}(:,length_Lsec(iii9,1):end) + (t_Vres*iii9); % add resin thickness so it moves towards x-direction
            y_left_parR{iii9,1} = y_left_par{end-iii9,1}(:,length_Lsec(iii9,1):end);
            
            % Values of 1st section x-axis
            xLsec{1,1}    = repmat(x_left_par{end,1}(:,:),length(Tot_layup),1);
            yLsec{1,1}    = repmat(y_left_par{end,1}(:,:),length(Tot_layup),1);
            
            % Store all other left sections
            xLsec{1,iii9+1}     = repmat(x_left_parR{iii9,1},length(Tot_layup),1);  
            yLsec{1,iii9+1}     = repmat(y_left_parR{iii9,1},length(Tot_layup),1); 
            
            % Locate where parent laminate starts from the x_rep
            idx_Rsec(iii9,1) = find(x_rep{1+iii9,:}(:,:) == idx_Rsec_temp(iii9,1));
            
            % For continuity of repair sections and right adhesive, take difference value of width in the repair and parent. This value will be added to right sections
            diff_xrepA = x_rep{1+1,:}(:,idx_Rsec(1,1)+1) - x_rep{1+1,:}(:,idx_Rsec(1,1)); % width difference in repair (it doesnt matter which layer we take)
            diff_xrepB = (x_rep{1+1,:}(:,idx_Rsec(1,1))) - xRrep_par(1,:); % width difference in repair and parent
            diff_xrep  = diff_xrepA - diff_xrepB;
            
            % Bottom repair layers used as reference for right resin and parent laminas; NOTE: diff_xrep is added to move all the values
            xRsec{1,1} = repmat(x_rep{1,1} + (t_Vres*length(L_remrow)),length(Tot_layup),1);
            yRsec{1,1} = repmat(y_rep{1,1},length(Tot_layup),1);
            
            % Store all other right sections, ensure same to add value of the move of th right sections after resin addition 
            xRsec{1,iii9+1}  = repmat((x_rep{1+iii9,:}(:,idx_Rsec(iii9,1):end))+ (t_Vres*iii9) + diff_xrep ,length(Tot_layup),1);
            yRsec{1,iii9+1}  = repmat(y_rep{1+iii9,:}(:,idx_Rsec(iii9,1):end),length(Tot_layup),1);    
        end
        
        % Ensure to add values to resins as well
        for ii9 = 1:length(L_remrow)
            xRrep1(ii9,1) = xres_right(ii9,1) + diff_xrep;
            yRrep1(ii9,1) = yres_right(ii9,1);
            
            % Build the x-axis which contains all resin
            xad_right{1,ii9}  = repmat(xRrep1(ii9,1),length(Tot_layup),1);
            yad_right{1,ii9}  = repmat(yRrep1(ii9,1),length(Tot_layup),1);
        end

        % Take last section which is isolated since we can just take the
        % end values of end right parent laminate. Ensure to add resin and
        % width difference.
        xRsecEnd = repmat(x_right_par{end,1}(:,:)+ (t_Vres*length(L_remrow) + diff_xrep),length(Tot_layup),1);
        yRsecEnd = repmat(y_right_par{end,1}(:,:),length(Tot_layup),1);
        
        % Put all cell into loop to concatetate them
        XLsec      = []; YLsec      = [];
        Xad_left   = []; Yad_left   = [];
        Xad_right  = []; Yad_right  = [];
        XRsec      = []; YRsec      = [];
        for i10  = 1:length(L_remrow)
            XLsec      = [ XLsec xLsec{1,i10} ]; 
            Xad_left   = [ Xad_left xad_left{1,i10} ]; 
            Xad_right  = [ Xad_right xad_right{1,i10} ]; 
            XRsec      = [ XRsec xRsec{1,i10}];
            
            YLsec      = [ YLsec yLsec{1,i10} ]; 
            Yad_left   = [ Yad_left yad_left{1,i10} ]; 
            Yad_right  = [ Yad_right yad_right{1,i10} ]; 
            YRsec      = [ YRsec yRsec{1,i10}];
        end

end