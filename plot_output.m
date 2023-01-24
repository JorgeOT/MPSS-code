%% Plot the output stress and strain from ABAQUS instead of the GUI
close all
namefile = 'Case_A';
[RF1,RF2,RM3,U1] = Read_ODB_outputs_node_function(namefile);
Totallengthofcompositemm = 615.1578990220405;
TotalThicknessmm = 6.9596;
plot_function(U1,Totallengthofcompositemm,RF1,TotalThicknessmm);

% % Create a table with the data and variable names
% T1 = table(RF1', U1','VariableNames', { '     RF1     ', '     U1     '} );
% T2 = table(Totallengthofcompositemm,TotalThicknessmm,E,'VariableNames',{'     Length     ','     Thickness   ','     E     '});
% % Write data to text file
% writetable(T1, 'RF1_U1.txt');
% writetable(T2, 'Length_Thickness_E.txt');
% 
