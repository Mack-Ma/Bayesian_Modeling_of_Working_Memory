%% SampleScript 4: Custom Model (Von Mises Distribution)
%
% Fit custom model (von mises distribution in this case)
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
ModelSpace={'VonMisesPDF'}; % Names of the custom models should perfectly match the function names
Nmodel=length(ModelSpace);
% Load toolbox
% addpath('Your BMW Toolbox Path')
BMW('Silent');
% Load Data
addpath('Your Data Path')
load('Your Data.mat') % File name should be revised accordingly
Nsubj=length(Data); % Variable name (Data) should be revised accordingly

%% Model Definition, Configuration & Estimation
MA_All=cell(1,Nmodel); % Pre-Allocation
for i=1:Nmodel
    fprintf('\n%s\n',ModelSpace{i})
    Config_MA.Data=Data;
    % Specify model
    Config_MA.Model.Model=ModelSpace{i};
    % Model Variants
    Config_MA.Model.Variants.Bias=0; % Deprecate bias here (bias=0)
    % Model definition
    MA=ModelDefinition_BMW(Config_MA);
    Config_MA.Criteria.Default=1;
%     Config_MA.Criteria.LLH=1; % Assign 1 to calculate log likelihood
%     Config_MA.Criteria.AIC=1; % Assign 1 to calculate AIC
%     Config_MA.Criteria.AICc=1; % Assign 1 to calculate AICc
%     Config_MA.Criteria.BIC=1; % Assign 1 to calculate BIC
%     Config_MA.Criteria.DIC=1; % Assign 1 to calculate DIC
%     Config_MA.Criteria.WAIC=1; % Assign 1 to calculate WAIC
%     Config_MA.Criteria.LME=1; % % Assign 1 to calculate log model evidence
    Config_MA.Constraints.start=10; % Initial guess
%     Config_MA.Constraint.start=[10, 0]; % When we use bias
    Config_MA.Constraints.lb=0;
%     Config_MA.Constraint.lb=[0, -90]; % When we use bias
    Config_MA.Constraints.ub=700;
%     Config_MA.Constraint.ub=[700, 90]; % When we use bias
    Config_MA.FitOptions.Algorithm='fmincon: interior-point'; % Optimization algorithm
    Config_MA.FitOptions.MaxIter=3000; % Max # of iteration
    Config_MA.FitOptions.Display='iter'; % Display mode
    Config_MA.InputFile=MA;
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
% Remove paths
rmpath('Your BMW Toolbox Path')
rmpath('Your Data Path')
