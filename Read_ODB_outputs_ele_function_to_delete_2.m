%function [S11,S22,S33,S12,LE11,LE22,LE33,LE12]=Read_ODB_outputs_ele_function(namefile,n)
clear all
%close all
namefile=('CASE_A');
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
    Flag = 0;
    tlinebefore = 'text';
    if(regexpi(tline, 'E L E M E N T   O U T P U T')>0)
        while(~contains(tline,'ERESIN')==1)
            tline = fgetl(fidd);
        end
    if (regexpi(tline, 'ERESIN')>0) %For elements, replace 'N O D E   O U T P U T'  by 'E L E M E N T   O U T P U T'
        
        numSteps = numSteps + 1;
%         for(re=1:12)
         tline = fgetl(fidd);
%         end
        i = i+1;
        %k=0;
%         while(isempty(str2num(tline)))
% %             if data_f(2) == 2 ||data_f(2) ==4
%                      tline = fgetl(fidd);
%         end
         while(isempty(str2num(tline)))
            if (contains(tline,'OR')==1 && contains(tline,'FOR')==0)
                table = readtable
                tline = replace(tline,'OR','');
                while(~isempty(str2num(tline)))
                    Flag = 1;
                    j=j+1;
                    data_f = sscanf(tline, '%d %e %e %e %e %e %e %e %e %e %e %e' , [1,10]);
                    data_f2(j,:)=data_f;
                    %             if (debug)
                    %                 fprintf(debug_file, '%d\t%e\n', data_f(1), data_f(2), data_f(3));
                    %             end
                    element_number=data_f(2);
                    i = i+1;
                    S11(j)=data_f(3);
                    S22(j)=data_f(4);
                    S33(j)=data_f(5);
                    S12(j)=data_f(6);
                    LE11(j)=data_f(7);
                    LE22(j)=data_f(8);
                    LE33(j)=data_f(9);
                    LE12(j)=data_f(10);
                    tlinebefore = tline;
                    tline = fgetl(fidd);
                    
                end
            end
            if (isempty(str2num(tline))==1)&&Flag==0
                tline = fgetl(fidd);
            end
            if contains(tline,'OR: *ORIENTATION USED FOR THIS ELEMENT')==1
                tline=fgetl(fidd);
            end
%             if isempty(str2num(tline)) && ~isempty(str2num(tlinebefore)) 
%                 tline = fgetl(fidd);
%             end
        end
            i = i+1;
    end
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
S11_output = S11;
S22_output = S22;
S33_output = S33;
S12_output = S12;
LE11_output = LE11;
LE22_output = LE22;
LE33_output = LE33;
LE12_output = LE12;
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
%end
