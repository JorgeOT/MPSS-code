function plot_function(U1,Totallengthofcompositemm,RF1,TotalThicknessmm)
figure();
%% For stress vs strain uncomment
%plot(U1/(Totallengthofcompositemm/1000),RF1/(TotalThicknessmm/1000));
% title('Stress strain');
% xlabel('Strain \epsilon');
% ylabel('Stress \sigma [N/m^2]');
%% For stress vs displacment uncomment
% plot(U1,RF1/(TotalThicknessmm/1000));
% title('Stress strain');
% xlabel('Displacement [m]');
% ylabel('Stress \sigma [N/m^2]');
%% For Load vs Strain uncomment
plot(U1*1000,RF1);
title('Stress strain');
xlabel('Displacement [mm]');
ylabel('Load [N]');

end