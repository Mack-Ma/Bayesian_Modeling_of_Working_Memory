%% SampleScript 4: Custom Model (Von Mises Distribution)
%
% Fit custom model(s) (von mises distribution in this case)
% ------------
% Sample script for BMW toolbox
% For detailed instruction of this script, please refer to the manual (type BMW('Manual')). 
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 12/8/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

%% Prologue
close all
clear variables
ModelSpace={'VonMisesPDF'}; % Names of the custom models should PERFECTLY match the function names
Nmodel=length(ModelSpace);
% Load toolbox & data path
addpath('Your toolbox path')
addpath('Your data path')
BMW('AddPath');

%% Model Configuration & Estimation
MA_All=cell(1,Nmodel); % Pre-Allocation
for i=1:Nmodel
    fprintf('\n%s\n',ModelSpace{i})
    Config_MA.Data='Your data directory'; % e.g. E:/matlab/Data_BMW.mat
    % Specify model
    Config_MA.Model.Model=ModelSpace{i};
    % Model Variants
    %     Config_MA.Model.Variants={'Bias'}; % Deprecate bias here (bias=0)
    % Specity mode
    Config_MA.FitOptions.Method='MLE';
    
    % The initial values and the constraints of the custom models need to be manually assigned
    Config_MA.Constraints.start=10; % Initial guess
    %     Config_MA.Constraint.start=[10, 0]; % When we use the location parameter (Bias)
    Config_MA.Constraints.lb=0;
    %     Config_MA.Constraint.lb=[0, -90]; % When we use the location parameter (Bias)
    Config_MA.Constraints.ub=700;
    %     Config_MA.Constraint.ub=[700, 90]; % When we use the location parameter (Bias)
    
    Config_MA.FitOptions.Algorithm='GA'; % Optimization algorithm
    % Configuration
    MA=Configuration_BMW(Config_MA);
    % Estimation
    MA=ModelFit_BMW(MA);
    MA_All{i}=MA;
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you‘re gonna save the results
save('Your Filename') % Save all variables
% Remove toolbox & data path
rmpath('Your toolbox path')
rmpath('Your data path')
