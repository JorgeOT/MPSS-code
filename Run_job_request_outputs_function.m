%function [RF1,RF2] = Run_job_request_outputs_function(namefile)%,read_output)

function [RF1,RF2,RM3,U1] = Run_job_request_outputs_function(namefile,N_cores)
% close all
S = mfilename('fullpath');
f = filesep;
ind=strfind(S,f);
S1=S(1:ind(end)-1);
cd(S1)
% namefile = 'Case_A';
% N_cores = 8;
%above sets the path
namefileodb = '.odb';
delete(strcat(namefile,namefileodb));
namefilelck = '.lck';
delete(strcat(namefile,namefilelck));
namefileinp = '.inp';
pause(2) % can this pause stop the job from getting stuck?
stringsystem1 = strcat('abaqus job=',namefile);
stringsystem2 = strcat(' cpus=',num2str(cast(N_cores,"int8")));
stringsystem3 = strcat(stringsystem1,' inp=',namefile,namefileinp,stringsystem2,' interactive');
%prueba = stringsystem2
system(stringsystem3); 



pause(2)
while exist(strcat(namefile,namefilelck),'file')==2
    pause(0.1)
end

while exist(strcat(namefile,namefileodb),'file')==0
    pause(0.1)
end

[RF1,RF2,RM3,U1]=Read_ODB_outputs_node_function(namefile);

end


%figure (1)
%plot(dis,force)