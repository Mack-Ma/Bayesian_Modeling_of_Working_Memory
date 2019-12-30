%% Sample Script 1: Model Fitting, Assessment & Comparison
%
% Assess and compare all candidate models of continuous recall
% ------------
% Sample script for BMW toolbox
% For detailed instruction of this script, please refer to the manual (type BMW('Manual')). 
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 10/31/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

%% Prologue
close all
clear variables
ModelSpace={'Item Limit', 'Standard Mixture', 'Slots-plus-Averaging', 'Equal Precision', 'Variable Precision', 'Variable Precision with Capacity',...
    'Categorical Slots-plus-Averaging', 'Categorical Variable Precision', 'Categorical Variable Precision with Capacity'}; % Full model space
Nmodel=length(ModelSpace);
% Load toolbox
addpath('Your BMW Toolbox Path')
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
    % Set Variants
    Config_MA.Model.Variants.Swap=1; % Swap rate
    Config_MA.Model.Variants.PrecF=0; % Fluctuation of precision
    Config_MA.Model.Variants.BiasF=0; % Fluctuation of bias
    Config_MA.Model.Variants.Bias=0; % Bias
    % Model definition
    MA=ModelDefinition_BMW(Config_MA);
%     Config_MA.Criteria.LLH=1; % Assign 1 to calculate log likelihood
%     Config_MA.Criteria.AIC=1; % Assign 1 to calculate AIC
%     Config_MA.Criteria.AICc=1; % Assign 1 to calculate AICc
%     Config_MA.Criteria.BIC=1; % Assign 1 to calculate BIC
%     Config_MA.Criteria.DIC=1; % Assign 1 to calculate DIC
%     Config_MA.Criteria.WAIC=1; % Assign 1 to calculate WAIC
%     Config_MA.Criteria.LME=1; % % Assign 1 to calculate log model evidence
    % All default
    Config_MA.Criteria.Default=1;  
    Config_MA.Constraints.Default=1;
    Config_MA.FitOptions.Default=1;
%     Config_MA.FitOptions.Algorithm='SA'; % Change optimization algorithm (Default: 'GA')
%     Config_MA.FitOptions.MaxIter=5000; % Change the max # of iteration (Default: 3000)
%     Config_MA.FitOptions.Display='off'; % Change the display mode (Default: 'iter')
    Config_MA.InputFile=MA;
    % Configuration
    MA=Configuration_BMW(Config_MA);
    % Estimation
    MA=ModelFit_BMW(MA);
    MA_All{i}=MA;
end

%% Model Comparison
Config_MC.MS=MA_All; % Define model space
Config_MC.Criterion='LLH';
MC_All.LLH=ModelComparison_BMW(Config_MC);
Config_MC.Criterion='AIC';
MC_All.AIC=ModelComparison_BMW(Config_MC);
Config_MC.Criterion='AICc';
MC_All.AICc=ModelComparison_BMW(Config_MC);
Config_MC.Criterion='BIC';
MC_All.BIC=ModelComparison_BMW(Config_MC);
% Config_MC.Criterion='LME'; % Log model evidence
% MC_All.LME=ModelComparison_BMW(Config_MC); % Make comparison based on LME

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
rmpath('Your Data Path')
