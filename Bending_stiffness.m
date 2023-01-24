function [D_rot] = Bending_stiffness(theta_bending1,RM3,length_composite)
% input from ABAQUS
            %% Bending stiffness from ABAQUS
            % Radious calculation
             theta_bending = 90-theta_bending1*180/pi;% in degrees
             centre_angle = 180-2*theta_bending;% in degrees
            % Sin theorem for triangles to calculate the radious
             R = length_composite*sin(theta_bending*pi/180)/(sin(centre_angle*pi/180));% (m)
            % Using rotational angle on 2D panel edges
            % Procedure : 1.) Take two points on each edges to construct two
            % lines to be intersected to find radius R, in order to
            % calculate curvature K = 1/R
             K_rot = 1/R; % [1/m]
             D_rot = RM3/ K_rot ; % [N.m^2]
end