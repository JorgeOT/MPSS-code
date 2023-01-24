%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: Trace all resins, parent laminate, repair laminate etc
% needed for ABAQUS input file generation.

function [Xres_el, Yres_el, rep_lay, parTrR, parTrL, hty_lay] = I_trace( X_res_Tot, Y_res_Tot,Xad_left,Xad_right, Nber_defect,Nber_nodefect,flag)

        Xad_left  = Xad_left - 1;
        Xad_right = Xad_right - 1;
        %% Detect layers that has defect
        id_def   = length(Y_res_Tot(:,1)) - Nber_defect + 1;
        xlay_def = X_res_Tot(id_def:end,:);
        ylay_def = Y_res_Tot(id_def:end,:);
        
        %% Trace resin in between parent and resins
        % Resin locations
        Xloc_adL = flip(Xad_left(1,:)); % bottom to top; % Note that Xad_left is matrix that is why we take only one row
        Xloc_adR = Xad_right(1,:); % Note that Xad_right is matrix that is why we take only one row
        
        % Take x and y values by considering different layers
        if (flag == 0) || (flag == 2)
            for i12 = 1:Nber_defect
                % Find idx location of resin in xlay_def on left
                xloc_adL(1,i12) = find(ismembertol(xlay_def(1,:)',Xloc_adL(1,i12),'ByRows',true));
                Yloc_adL(1,i12) = ylay_def(i12,xloc_adL(1,i12)); % corresponidng y values
                
                % Find idx location of resin in xlay_def on right
                xloc_adR(1,i12) = find(ismembertol(xlay_def(1,:)',Xloc_adR(1,i12),'ByRows',true));
                Yloc_adR(1,i12) = ylay_def(i12,xloc_adR(1,i12));  % corresponidng y values
                
                % Take neighbour (point before) x and y values to close element of resin
                xnbh_adL(1,i12) = xloc_adL(1,i12)-1; % take previous point
                xnbh_adR(1,i12) = xloc_adR(1,i12)-1; % take previous point

                % Store neighbour point
                Xloc_adL0(1,i12) = xlay_def(1, xnbh_adL(1,i12));
                Xloc_adR0(1,i12) = xlay_def(1, xnbh_adR(1,i12));
                Yloc_adL0(1,i12) = ylay_def(i12,xnbh_adL(1,i12));
                Yloc_adR0(1,i12) = ylay_def(i12,xnbh_adR(1,i12));  
            end
            
            %% Trace the repair layers with new x-axis
            % Find loc indx of the in and end of repair layers at X and corres.Y
            for i = 1:Nber_defect
                idx_repIn(i,1)   = find(ismembertol(xlay_def(1,:)', Xloc_adL(1,i),'ByRows',true));           
                idx_repIEnd(i,1) = find(ismembertol(xlay_def(1,:)', Xloc_adR0(1,i),'ByRows',true)) - 1; % -1 to not include resins on right
                xrep_lay{i,:}    = xlay_def(i,idx_repIn(i,1):idx_repIEnd(i,1))';
                yrep_lay{i,:}    = ylay_def(i,idx_repIn(i,1):idx_repIEnd(i,1))';
                rep_lay{i,1}     = [xrep_lay{i,:} yrep_lay{i,:}];
            end
            
            %% Trace the parent layers with new x-axis from Nber_defect
            for i = 1:Nber_defect
                idx_parLend(i,1) = find(ismembertol(xlay_def(1,:)', Xloc_adL0(1,i),'ByRows',true));
                idx_parRin(i,1)  = find(ismembertol(xlay_def(1,:)', Xloc_adR(1,i),'ByRows',true));
                x_parTrL{i,:} = xlay_def(i,1:idx_parLend(i,1)-1)';% -1 to not include resins on left
                y_parTrL{i,:} = ylay_def(i,1:idx_parLend(i,1)-1)';
                x_parTrR{i,:} = xlay_def(i,idx_parRin(i,1):end)';
                y_parTrR{i,:} = ylay_def(i,idx_parRin(i,1):end)';
                % store in pairs
                parTrR{i,1}   = [x_parTrR{i,:} y_parTrR{i,:}];
                parTrL{i,1}   = [x_parTrL{i,:} y_parTrL{i,:}];
            end
            
            else if (flag == 1) || (flag == 3)
                    for i12 = 1:Nber_defect
                        % Find idx location of resin in xlay_def on left
                        xloc_adL(1,i12) = find(round(xlay_def(1,:), 4) == round(Xloc_adL(1,i12), 4));
                        Yloc_adL(1,i12) = ylay_def(i12,xloc_adL(1,i12)); % corresponidng y values
                        
                        % Find idx location of resin in xlay_def on right
                        xloc_adR(1,i12) = find(round(Xloc_adR(1,i12), 4) == round(xlay_def(1,:), 4));
                        Yloc_adR(1,i12) = ylay_def(i12,xloc_adR(1,i12));  % corresponidng y values
                        
                        % Take neighbour (point before) x and y values to close element of resin
                        xnbh_adL(1,i12) = xloc_adL(1,i12)-1; % take previous point
                        xnbh_adR(1,i12) = xloc_adR(1,i12)-1; % take previous point
            
                        % Store nneighbour point
                        Xloc_adL0(1,i12) = xlay_def(1, xnbh_adL(1,i12));
                        Xloc_adR0(1,i12) = xlay_def(1, xnbh_adR(1,i12));
                        Yloc_adL0(1,i12) = ylay_def(i12,xnbh_adL(1,i12));
                        Yloc_adR0(1,i12) = ylay_def(i12,xnbh_adR(1,i12));
                    end
                    
                    %% Trace the repair layers with new x-axis

                    % find loc indx of the in and end of repair layers at X and corres.Y
                    for i = 1:Nber_defect
                        idx_repIn(i,1)   = find(round(xlay_def(i,:), 4) == round(Xloc_adL(1,i), 4));
                        idx_repIEnd(i,1) = find(round(xlay_def(i,:), 4) == round(Xloc_adR0(1,i), 4)) - 1; % -1 to not include resins on right
                        xrep_lay{i,:}    = xlay_def(i,idx_repIn(i,1):idx_repIEnd(i,1))';
                        yrep_lay{i,:}    = ylay_def(i,idx_repIn(i,1):idx_repIEnd(i,1))';
                        rep_lay{i,1}     = [xrep_lay{i,:} yrep_lay{i,:}];
                    end
                    
                    %% Trace the parent layers with new x-axis from Nber_defect
                    for i = 1:Nber_defect
                        idx_parLend(i,1) = find(round(xlay_def(i,:), 4) == round(Xloc_adL0(1,i), 4));
                        idx_parRin(i,1)  = find(round(xlay_def(i,:), 4) == round(Xloc_adR(1,i), 4));
                        x_parTrL{i,:} = xlay_def(i,1:idx_parLend(i,1)-1)';% -1 to not include resins on left
                        y_parTrL{i,:} = ylay_def(i,1:idx_parLend(i,1)-1)';
                        x_parTrR{i,:} = xlay_def(i,idx_parRin(i,1):end)';
                        y_parTrR{i,:} = ylay_def(i,idx_parRin(i,1):end)';
                        % store in pairs
                        parTrR{i,1}   = [x_parTrR{i,:} y_parTrR{i,:}];
                        parTrL{i,1}   = [x_parTrL{i,:} y_parTrL{i,:}];
                    end
                    
                    end
        end

        % Put all resin points into one vector
        xres1 = [ flip(Xloc_adL0) Xloc_adR0 ];
        yres1 = [ flip(Yloc_adL0) Yloc_adR0 ];
        xres2 = [ flip(Xloc_adL)  Xloc_adR  ];
        yres2 = [ flip(Yloc_adL)  Yloc_adR  ];
        
        % Combine all
        Xres_el = [ xres1 xres2]';
        Yres_el = [ yres1 yres2]';
                
       %% Trace the healthy layers with new x-axis from Nber_nodefect
        Hty_lay = [];
        for i = 1:Nber_nodefect
            xlay_nodef(:,i) = X_res_Tot(i,:);
            ylay_nodef(:,i) = Y_res_Tot(i,:);
            hty_lay{:,i}    = [xlay_nodef(:,1) ylay_nodef(:,i)];
        end

end