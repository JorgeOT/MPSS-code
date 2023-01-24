%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

% Description: Take the half matrix of defected layers in order to have to cuts of panel (left and right
% cuts). Find the mid-value of defect width and half length of defect, so
% we know where to start the scarfing (i.e. removal of layers according to material, 
% i.e. 100 mm for UD and 50 mm for BIAX). The 'cut' matrices will contain
% two parts: 'mid' and 'pos' matrices. 'Pos' matrix takes the matrix where
% removal of scarf length will be performed. 'mid' matrix takes the other part of 'cut' matrix which 
% contains defect according to half-length of defect

function [xrep_left_aseb, yrep_left_aseb, xrep_right_aseb, yrep_right_aseb,x_right_par, y_right_par, x_left_par, y_left_par, L_remrow] = D_removedefect(layup, flag,x_aseb, y_aseb, x_defect, Nber_nodefect, layup_defect, ud, biax, xud_olap, xbiax_olap)
                                            
        % Take x and y from new data where bottom layer is at y = 0
        x_stores = x_aseb(Nber_nodefect+1:end,:);
        y_stores = y_aseb(Nber_nodefect+1:end,:);
        
        % Find the center distance at x direction and cut section for right part
        loc_center = length(x_stores)/2;
        xcut_right  = x_stores(:, loc_center:end);
        ycut_right  = y_stores(:, loc_center:end);                          % Bottom to top direction

        cnt_right = 1;                                                      % counter
        
        % Removal of defect on right
        mid_defect = x_defect/2;  % Half of defect width
        
        % Width of defect from center from first lamina at the bottom layer
        for i4 = 1: loc_center
            % Subtract next element value to previous value to get difference
            xdiff_right = round(xcut_right(1, i4+1) - xcut_right(1, 1),2);
            cnt_right     = i4 + 1;                                         % Index of mid-witdth of defect
            
            % If the difference is higher than the half-length of defect them we take the matrix where we initiate the scarfing
            if xdiff_right > mid_defect
                % 'Pos' matrix takes the matrix where removal of scarf length will be performed
                xpos_right = xcut_right(:, cnt_right:end);
                ypos_right = ycut_right(:, cnt_right:end); 

                % 'mid' matrix takes the other part of 'cut' matrix which contains defect according to half-length of defect
                xmid_right = xcut_right(:, 1:cnt_right);
                ymid_right = ycut_right(:, 1:cnt_right);
                break;
            else
                continue;
            end
        end
        
        % Take layup with direction top to bottom
        layup_rem  = flip(layup_defect);  
        
        % Take layup size to store location of scarf length removal according to material of the lamina (UD or Biax)
        layup_Rem  = zeros(size(layup_rem));

        if(flag == 0)
            % Take index location of UD and Biax material
            idUd_rem   = find(layup_rem == ud);
            idBiax_rem = find(layup_rem == biax);       

            % Apply length of scarfing to remove, put into a vector containing scarf length info
            layup_Rem(idUd_rem,:)   = xud_olap; 
            layup_Rem(idBiax_rem,:) = xbiax_olap; 
            
            % Create summed vector to sum-up previous scarf length to next one
            Layup_Rem = sum(layup_Rem) * ones(size(layup_Rem));

            % From the summed vector, subtract the vector containing scarf length info
            for i5a = 1:length(layup_Rem)-1
                Layup_Rem(end-i5a,:) = Layup_Rem(end-i5a+1 ,:) - layup_Rem(end-i5a+1);
            end

            % According to layer material type, then add the corresponding
            % scarf length by material in the xpos_right matrix according
            % to layer number
            for i6a = 1 : length(Layup_Rem)
                countidUd = i6a;
                if any(idUd_rem(:) == countidUd)
                    udlay_rem(countidUd,:) = xpos_right(countidUd,1) + Layup_Rem(countidUd); 
                        else if any(idBiax_rem(:) == countidUd)
                                biaxlay_rem(countidUd,:) = xpos_right(countidUd,1) + Layup_Rem(countidUd); 
                            end  
                end
            end      

            % Taking non zero indeces where UD and Biax layers are located
            UDindex = exist ('udlay_rem','var');
            BIAXindex = exist ('biaxlay_rem','var');
            if UDindex == 1
            udrow    = find(udlay_rem~=0);
            end
            if BIAXindex == 1
            biaxrow  = find(biaxlay_rem~=0);
            end
            
            if (UDindex + BIAXindex) == 2
            L_remrow = sort([ udlay_rem(udrow); biaxlay_rem(biaxrow) ]);
            elseif (UDindex) == 1
                L_remrow = sort([ udlay_rem(udrow)]);
            elseif (BIAXindex) == 1
                L_remrow = sort([biaxlay_rem(biaxrow) ]);
            end
            
            %% Skip the other flags and proceed to line 128
            
            else if (flag == 1)
                % Apply length of scarfing to remove, put into a vector
                if layup(1,1) == biax
                    idBiax_rem = find(layup_rem == biax);
                    layup_Rem(idBiax_rem,:) = xbiax_olap;
                    else if layup(1,1) == ud
                        idUd_rem   = find(layup_rem == ud);
                        layup_Rem(idUd_rem,:)   = xud_olap; % [100 100 100 100 100]
                        end
                end
                Layup_Rem = sum(layup_Rem) * ones(size(layup_Rem)); 

                % Create a vector of length to remove by subtracting next layer, starting from total
                for i5b = 1:length(layup_Rem)-1
                    Layup_Rem(end-i5b,:) = Layup_Rem(end-i5b+1 ,:) - layup_Rem(end-i5b+1); 
                end

                    % According to layer type, then add the scarfing length
                    for i6b = 1 : length(Layup_Rem)
                        countidUd = i6b;
                        if layup(1,1) == ud
                            udlay_rem(countidUd,:) = xpos_right(countidUd,1) + Layup_Rem(countidUd);

                            % Taking non zero indeces where UD layers are located from vector
                            udrow = find(udlay_rem~=0);
                            L_remrow = sort(udlay_rem(udrow));

                                else if layup(1,1) == biax
                                        biaxlay_rem(countidUd,:) = xpos_right(countidUd,1) + Layup_Rem(countidUd);

                                        % Taking non zero indeces where biax layers are located from vector
                                        biaxrow  = find(biaxlay_rem~=0);
                                        L_remrow = sort(biaxlay_rem(biaxrow));
                                    end
                        end
                    end
                end
        end

        % Store all parent laminate on the right
        for i7 = 1 : length(L_remrow)
            % Take location of L_remrom on the 'pos' matrix
            % This 'if' statements ensures to find index location for any
            % case
            tol = 1;% tolerance for the find and the ismember command
            if ismember(round(L_remrow), round(xpos_right(1,:)))
                    id_rem     = find(round(xpos_right(i7,:)) == round(L_remrow(i7,:)));
            else if ismember(fix(L_remrow), fix(xpos_right(1,:)))
                    id_rem     = find(fix(xpos_right(i7,:)) == fix(L_remrow(i7,:)));
                else if ismember(L_remrow, xpos_right(1,:))
                    id_rem     = find(xpos_right(i7,:) == L_remrow(i7,:));
                    else if ismembertol(L_remrow, xpos_right(1,:),tol)
                    id_rem     = find(abs(xpos_right(i7,:)-L_remrow(i7,:))<tol);
                        end
                            %id_rem     = find(xpos_right(i7,:) == L_remrow(i7,:));
                    end
                end
            end
            
            x_right_par{i7,1} = xpos_right(i7,id_rem:end);
            y_right_par{i7,1} = ypos_right(i7,id_rem:end);

            % Taking length from mid axis to overlap on ith layer
            Xdiff_right{i7,1} = xpos_right(i7,1:id_rem); 
            Ydiff_right{i7,1} = ypos_right(i7,1:id_rem);      

            % Assemble repair layers of right
            xrep_right_aseb{i7,1} = [xmid_right(i7,:), Xdiff_right{i7,1}];
            yrep_right_aseb{i7,1} = [ymid_right(i7,:), Ydiff_right{i7,1}];
        end

        %% Removal of defect on Left
        % Since the size of the left size is known, we can just take it from right ensuring symmetry
         xcut_left = fliplr(x_stores(:,1:loc_center));
         ycut_left = fliplr(y_stores(:,1:loc_center));                        

            for i8 = 1: length(L_remrow)
                % Parent laminate at right
                x_left_par{i8,1} = fliplr(xcut_left(i8,length(xrep_right_aseb{i8,1}):end));
                y_left_par{i8,1} = fliplr(ycut_left(i8,length(yrep_right_aseb{i8,1}):end));

                % Assemble defect layers to be removed for right
                xrep_left_aseb{i8,1} = fliplr(xcut_left(i8,(1:length(xrep_right_aseb{i8,1}))));
                yrep_left_aseb{i8,1} = fliplr(ycut_left(i8,(1:length(yrep_right_aseb{i8,1})))); 
            end

end