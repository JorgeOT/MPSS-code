% Xiao Chen, Wilfried Njomo-Wandji, Xing-Yuan Miao (2022), 
% A robust and automated method for geometric modelling of thick laminates with multiple and asymmetric ply wrinkles
% https://doi.org/10.1016/j.compstruct.2022.115319, Elsevier Ltd., Composites Structures

% Multi-Parameter Stepped Scarf (MPSS)
% Modified by: Aura Venessa Paguagan (venessa.paguagan@outlook.com)

function [x_store, y_store] = B_peakfitting(flag, Wrinkle, Nber_defect, disc_density, centre_guess, sigma_guess, amplitude_guess, layup_defect, t_Hres)

%% Initialise Gaussian not constant parameters - a_ij, b_ij, c_ij 
numGaussians = zeros(1,Nber_defect);
meanResidual = zeros(1,Nber_defect);
centres = cell(1,Nber_defect);
sigmas = cell(1,Nber_defect);
amplitudes = cell(1,Nber_defect);
tActual = cell(1,Nber_defect);
x_store = zeros(Nber_defect,disc_density);                                  % Store discretized x
y_store = zeros(Nber_defect,disc_density);                                  % Store discretized y

xy = Wrinkle{1,1};                                                          % xy access each double at each cell
unique(xy.','rows').';                                                      % Sort the data

%Inspected Data points into two columns
X = xy(:,1);
Y = xy(:,2);
[x,sortIndex] = sort(X);                                            % taking X and index row number
y = Y(sortIndex);

while 1
[pks,idxs,widths,proms] = findpeaks(y,'MinPeakProminence',0.1,'Annotate','extents'); % Find indices of the peaks. 'idxs' is peak location.

    % Let N_i^Gauss = N_i^peak (pks):
    if numel(pks)< 3
       numGaussians = numel(pks)+5;
       centres{1,1} = randi([1,centre_guess],1,numGaussians);
       sigmas{1,1} = randi([1,sigma_guess],1,numGaussians);
       amplitudes{1,1} = randi([1,amplitude_guess],1,numGaussians);              
    else if numel(pks)< 9
            numGaussians = numel(pks)+7;
            centres{1,1} = randi([1,centre_guess],1,numGaussians);
            sigmas{1,1} = randi([1,sigma_guess],1,numGaussians);
            amplitudes{1,1} = randi([1,amplitude_guess],1,numGaussians);
        else if numel(pks)< 15
                numGaussians = numel(pks)+9;
                centres{1,1} = randi([1,centre_guess],1,numGaussians);
                sigmas{1,1} = randi([1,sigma_guess],1,numGaussians);
                amplitudes{1,1} = randi([1,amplitude_guess],1,numGaussians);
    else
        numGaussians = numel(pks);
        centres{1,1} = randi([1,centre_guess],1,numGaussians);
        sigmas{1,1} = randi([1,sigma_guess],1,numGaussians);
        amplitudes{1,1} = randi([1,amplitude_guess],1,numGaussians);                
            end
        end
    end

% Create a table gathering all the parameters
tActual{1,1} = table((1:numGaussians(1))',amplitudes{1,1}(:),...
    centres{1,1}(:),sigmas{1,1}(:),'VariableNames',{'Number','Amplitude','Mean','Width'});
% Sort parameters with order of increasing mean
tActual{1,1} = sortrows(tActual{1,1},3);
tActual{1,1}.Number = (1:numGaussians(1))'; % Unsort the first column of numbers: Gaussian

% Fit Gaussian peaks:
initialGuesses  = [tActual{1,1}.Mean(:),tActual{1,1}.Width(:)];
startingGuesses = reshape(initialGuesses',1,[]);

global c
x = reshape(x,1,[]);                                            % data points in length
y = reshape(y,1,[]);                                            % data points in height

% Perform an iterative fit using the FMINSEARCH
options = optimset('TolFun',1e-4,'TolX',1e-4,'MaxIter',10^8,'MaxFunEvals',10^8);

[parameter,fval,exitflag,output] = fminsearch(@(lambda)(fitgauss(lambda,x,y))...
    ,startingGuesses,options); 

% Calculate through Gaussian function
yhat = PlotComponentCurves(x,c,parameter);

% Relative and Residual error check
    %'if' the mean relative error ('meanResiduals') between y_i(actual y) and y_i^fit (estimated yhat) is smaller than eta error(0.05)
    % and the relative error ('peakResidual') between each local peak y_j^peak and y_j^fit,peak is smaller thea error 0.05
meanResidual = mean(abs(y-yhat)/y);

% Control ends error
Rend = abs(y(end)-yhat(end))/y(end);
Lend = abs(y(1)-yhat(1))/y(1);

% Compute the peak residuals between the atual y and estimated y
for i5a = 1:size(idxs)
    idx = idxs(i5a);
    peakResidual(1,i5a) = abs(y(idx)-yhat(idx))/y(idx);
end

% Satisfy both mean error and peak error criterion
if meanResidual<0.05 && all(peakResidual(1,:)<0.05) && Rend<0.05 && Lend<0.05
    break
end

end % 2nd 'while 1

% Gather all parameters for Gaussian function
estimatedMuSigma   = reshape(parameter,2,[])';
gaussianParameters = [c,estimatedMuSigma];
gaussianParameters = sortrows(gaussianParameters,2);

%% Discretize data
x_lin        = linspace(min(x),max(x),disc_density);
X_store(1,:) = x_lin;                                       
y_lin        = PlotComponentCurves(x_lin,c,parameter);         % Obtain corresponding discretized y_i^fit,disc
Y_store(1,:) = y_lin;
        
if(flag == 0) % Composite with Defect - if layup contains two different fibre material
    % Create matrix of defect
    x_Wrinkles = repmat(X_store, Nber_defect,1);
    y_Wrinkle  = repmat(Y_store, Nber_defect,1);
    
    % Create temporary layer which contains layup information
    temp_lays = layup_defect.* ones(length(layup_defect), length(y_Wrinkle)); 
    temp_lay = [zeros(1, length(y_Wrinkle)); temp_lays];
    
    for ia = 1 : Nber_defect
        % Add the previous value to the next one, until values are accumulated
        temp_lay(ia + 1,:) = temp_lay(ia + 1,:) +  temp_lay(ia ,:);
        
        % Bottom layers with defect
        y_Diff(ia,:) = y_Wrinkle(ia,:) - temp_lay(ia+1,:); 
    end
    y_Wrinkles = flip([ Y_store ; y_Diff(1:end-1,:)]); % top to bottom direction

else if (flag == 1) % Composite with Defect - if layup contains only one type of fibre material
        % Create matrix of defect
        x_Wrinkles = repmat(X_store, Nber_defect,1);
        y_Wrinkle  = repmat(Y_store, Nber_defect,1);
        
        % Create temporary layer which contains layup information
        temp_lays = layup_defect.* ones(length(layup_defect), length(y_Wrinkle)); 
        temp_lay = [zeros(1, length(y_Wrinkle)); temp_lays];
        for ib = 1 : Nber_defect
            temp_lay(ib + 1,:) = temp_lay(ib + 1,:) +  temp_lay(ib ,:);
            y_Diff(ib,:) = y_Wrinkle(ib,:) - temp_lay(ib+1,:); 
        end
        y_Wrinkles = flip([ Y_store ; y_Diff(1:end-1,:)]);
    end
end

x_store = x_Wrinkles;
y_store = y_Wrinkles;
y_store(end,:) = y_store(end,:) - t_Hres;

% Plot Figure 
% Figure size
x0     = 50;
y0     = 500;
width  = 2000;
height = 400;

figure()
plot(x,yhat,'.r','LineWidth',1,'MarkerSize',11);
xlabel('length / mm','Interpreter','latex','Fontsize',13);
ylabel('thickness / mm','Interpreter','latex','Fontsize',13);
set(gcf,'position',[x0,y0,width,height])
xlim(sort([x(1) x(end)]))

figure()
plot(x_store(1,:),y_store,'-k','LineWidth',1,'MarkerSize',11);
xlabel('length / mm','Interpreter','latex','Fontsize',13);
ylabel('thickness / mm','Interpreter','latex','Fontsize',13);
set(gcf,'position',[x0,y0,width,height])
xlim(sort([x(1) x(end)]))
        
end

function yhat = PlotComponentCurves(t,c,parameter) % yhat is an estimated y
try
    % Get the means and widths. 
    means = parameter(1:2:end);
    widths = parameter(2:2:end);
    yhat = zeros(1,length(t));
    numGaussians = length(c);
    legendStrings = cell(2,1);
    for k=1:numGaussians
        % Get each component curve
        thisEstimatedCurve = c(k).*gaussian(t,means(k),widths(k));
        yhat = yhat + thisEstimatedCurve;
        %legendStrings{k} = sprintf('Estimated Gaussian %d', k);
    end
end
end

%==========================================================================================================================================
function theError = fitgauss(lambda,t,y)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the % way, for educational purposes.
% Author: T. C. O'Haver, 2006

global c 
try
    A = zeros(length(t),round(length(lambda)/2));
    for j = 1:length(lambda)/2
       A(:,j) = gaussian(t,lambda(2*j-1),lambda(2*j))';
    end
    c = A\y';
    z = A*c;
    theError = norm(z - y');
end
end

%==========================================================================================================================================
function gauss = gaussian(x,peakPosition,width) % (a,b,c gaussian parameters)
% gaussian(x,position,width) = gaussian peak centred on position, half-width=width
% x may be scalar, vector, or matrix, position and width both scalar 
% Author: T. C. O'Haver, 1988
gauss = exp(-((x-peakPosition)./width).^2/2) ./ (sqrt(2*pi)*width);

end
