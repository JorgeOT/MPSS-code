function [S11,S22,S33,S12,LE11,LE22,LE33,LE12]=Read_ODB_outputs_ele_function2(namefile,elset)
%clear all
%close all
%namefile=('CASE_A');
debug = 0;
if (debug)
    debug_file = fopen('debug.dat', 'wt');
end
namefiledat = '.dat';
while exist(strcat(namefile,namefiledat),'file')==0
    pause(0.1)
end
file = strcat(namefile,namefiledat);
fidd = fopen(file);
i = 0;
pause(0.5)
numSteps = 0;
j=0;
data_f=zeros(1,4);
while ( ~feof(fidd) )
    
    tline = fgetl(fidd);
    
    i = i+1;

    if(regexpi(tline, 'E L E M E N T   O U T P U T')>0)
        while(~contains(tline,elset)==1)
            tline = fgetl(fidd);
        end
        while(~contains(tline,'MAXIMUM')==1)
            tline = fgetl(fidd);
        end
        if (regexpi(tline, 'MAXIMUM')>0) %For elements, replace 'N O D E   O U T P U T'  by 'E L E M E N T   O U T P U T'
                    tline = replace(tline,'MAXIMUM','');
                    data_f = sscanf(tline, '%e %e %e %e %e %e %e %e %e %e' , [1,8]);                  
        end
        while(~contains(tline,'MINIMUM')==1)
            tline = fgetl(fidd);
        end
        if (regexpi(tline, 'MINIMUM')>0) %For elements, replace 'N O D E   O U T P U T'  by 'E L E M E N T   O U T P U T'
                    tline = replace(tline,'MINIMUM','');
                    data_f2 = sscanf(tline, '%e %e %e %e %e %e %e %e %e %e' , [1,8]);                  
        end

        if (debug)
            fprintf(debug_file, '\n');
        end
    end
    
end

S11 = [data_f(1),data_f2(1)];
S22 = [data_f(2),data_f2(2)];
S33 = [data_f(3),data_f2(3)];
S12 = [data_f(4),data_f2(4)];
LE11 = [data_f(5),data_f2(5)];
LE22 = [data_f(6),data_f2(6)];
LE33 = [data_f(7),data_f2(7)];
LE12 = [data_f(8),data_f2(8)];
% dimensionS11 = size(S11);
% S11_output = S11(1,(dimensionS11(1,2)-(n*8-1)):end);
% S22_output = S22(1,(dimensionS11(1,2)-(n*8-1)):end);
% S33_output = S33(1,(dimensionS11(1,2)-(n*8-1)):end);
% S12_output = S12(1,(dimensionS11(1,2)-(n*8-1)):end);
% LE11_output = LE11(1,(dimensionS11(1,2)-(n*8-1)):end);
% LE22_output = LE22(1,(dimensionS11(1,2)-(n*8-1)):end);
% LE33_output = LE33(1,(dimensionS11(1,2)-(n*8-1)):end);
% LE12_output = LE12(1,(dimensionS11(1,2)-(n*8-1)):end);
if (debug)
    fclose(debug_file);
end
fclose(fidd);
fclose all
end
