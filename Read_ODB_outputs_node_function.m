function [RF1,RF2,RM3,U1]=Read_ODB_outputs_node_function(namefile)
%clear all
%close all

debug = 0;
if (debug)
    debug_file = fopen('debug.dat', 'wt');
end
% namefile = CASE_A;
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
while ( ~feof(fidd) )
    
    tline = fgetl(fidd);
    
    i = i+1;
    if (regexpi(tline, 'N O D E   O U T P U T')>0)  %For elements, replace 'N O D E   O U T P U T'  by 'E L E M E N T   O U T P U T'
        numSteps = numSteps + 1;
        tline = fgetl(fidd);
        i = i+1;
        
        while(isempty(str2num(tline)))
            tline = fgetl(fidd);
            i = i+1;
        end
        while(~isempty(str2num(tline)))
            j=j+1;
            data_f = sscanf(tline, '%d %e %e %e %e', [1,5]);
            if (debug)
                fprintf(debug_file, '%d\t%e\n', data_f(1), data_f(2), data_f(3), data_f(4), data_f(5));
            end
            node_number=data_f(1);
            tline = fgetl(fidd);
            i = i+1;
            RF1(j)=data_f(2);
            RF2(j)=data_f(3);
            RM3(j)=data_f(4);
            U1(j)=data_f(5);
            
        end
        if (debug)
            fprintf(debug_file, '\n');
        end
    end
end
if (debug)
    fclose(debug_file);
end
fclose(fidd);
fclose all;
end
