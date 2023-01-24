%function [D_sim,F_sim]=Read_ODB_outputs_ele_function(namefile)
clear all
%close all
debug = 0;
namefile = 'CASE_A';
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
    
    if (regexpi(tline, 'E L E M E N T   O U T P U T')>0)  %For elements, replace 'N O D E   O U T P U T'  by 'E L E M E N T   O U T P U T'
        
        numSteps = numSteps + 1;
        for(re=1:12)
        tline = fgetl(fidd);
        end
        i = i+1;
        k=0;
        while(isempty(str2num(tline)))
%             if data_f(2) == 2 ||data_f(2) ==4
%                     tline = fgetl(fidd);
%             end
            if (contains(tline,'OR')==1 && contains(tline,'FOR')==0)
                tline = replace(tline,'OR','');
                while(~isempty(str2num(tline)))
                    j=j+1;
                    data_f = sscanf(tline, '%d %e %e %e' , [1,4]);
                    data_f2(j,:)=data_f;
                    %             if (debug)
                    %                 fprintf(debug_file, '%d\t%e\n', data_f(1), data_f(2), data_f(3));
                    %             end
                    element_number=data_f(2);
                    tline = fgetl(fidd);
                    i = i+1;
                    D_sim(j)=data_f(3);
                    F_sim(j)=data_f(4);
                    
                    k=k+1;
                end
            end
            if k==32
                tline=fgetl(fidd);
            end
        end
            i = i+1;
        end
%         while(~isempty(str2num(tline)))
%             j=j+1;
%             data_f = sscanf(tline, '%d %e %s %e %e' , [1,5]);
% %             if (debug)
% %                 fprintf(debug_file, '%d\t%e\n', data_f(1), data_f(2), data_f(3));
% %             end
%             element_number=data_f(2);
%             tline = fgetl(fidd);
%             i = i+1;
%             %D_sim(j)=data_f(4);
%             %F_sim(j)=data_f(5);
% 
%         end
        if (debug)
            fprintf(debug_file, '\n');
        end
end
dimensionD = size(D_sim);
D_sim2 = D_sim(1,(dimensionD(1,2)-31):end);
dimensionF = size(F_sim);
F_sim2 = F_sim(1,(dimensionF(1,2)-31):end);
if (debug)
    fclose(debug_file);
end
fclose(fidd);
fclose all
%end
