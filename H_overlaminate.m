%% Multi-Parameter Stepped Scarf (MPSS)
% Author: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

function [X_res_Tot, Y_res_Tot, y_over, xy_over] = H_overlaminate(x_res_Tot, y_res_Tot, t_overmat, t_Hres, over_olap, rep_lay, xover_olap)

y_res_Tot(end,:) = y_res_Tot(end,:) + t_Hres;
y_over    = y_res_Tot(end,:) + t_overmat - t_Hres;

% If user wants to define overlapping length for overlaminate layer
if over_olap == 1
    top_in     = rep_lay{end,1}(1,1) - xover_olap;
    top_end    = rep_lay{end,1}(end,1) + xover_olap;
    tolin = 2;
    if ismember(round(top_in), round(x_res_Tot(1,:)))
        loc_topin  = find(round(top_in)== round(x_res_Tot(1,:)));
    else if ismember(fix(top_in), fix(x_res_Tot(1,:)))
            loc_topin  = find(fix(top_in)== fix(x_res_Tot(1,:)));
        else if ismember(top_in, x_res_Tot(1,:))
                loc_topin  = find(x_res_Tot(1,:) == top_in);
            else if ismembertol(top_in, x_res_Tot(1,:),tolin)
                loc_topin  = find(abs(x_res_Tot(1,:)-top_in)<tolin);
                end
            end
        end
    end
    tolend = 5;
    if ismember(round(top_end), round(x_res_Tot(1,:)))
        loc_topend  = find(round(top_end)== round(x_res_Tot(1,:)));
    else if ismember(fix(top_end), fix(x_res_Tot(1,:)))
            loc_topend  = find(fix(top_end)== fix(x_res_Tot(1,:)));
        else if ismember(top_end, x_res_Tot(1,:))
                loc_topend  = find(x_res_Tot(1,:) == top_end); 
            else if ismembertol(top_end, x_res_Tot(1,:),tolend)
                loc_topend  = find(abs(x_res_Tot(1,:)-top_end)<tolend);
                end
            end
        end
    end    
    x_over    = x_res_Tot(1,loc_topin:loc_topend);   
    
    % Left top layer to be the same without length change to get assembly
    X_res_Tot = [x_res_Tot; x_res_Tot(end,:)];
    Y_res_Tot = [y_res_Tot; y_over];
    
    % To trace top layer later on in ABaqus.m function
    xy_over   = [x_over' y_over(1,loc_topin:loc_topend)'];
    
    else
        X_res_Tot = [x_res_Tot; x_res_Tot(end,:)];
        Y_res_Tot = [y_res_Tot; y_over];
        xy_over = [x_res_Tot(end,:)' y_over'];
end
end