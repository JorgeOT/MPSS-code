function [E] = Youngs_modulus(U1,RF1,r,t)
% input from ABAQUS
            %% Young modulus E from ABAQUS results in GPa         
             E = (RF1*r)/(U1*t)/(10^6);
end