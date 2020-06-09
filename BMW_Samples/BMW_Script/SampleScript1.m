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
% Sun Yat-Sen University
% 10/31/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

%% Prologue
close all
clear variables
ModelSpace={'Item Limit', 'Standard Mixture', 'Slots-plus-Averaging', 'Equal Precision', 'Variable Precision', 'Variable Precision with Capacity',...
    'Category-Only','Category-Only (with Capacity)'}; % Full model space
Nmodel=length(ModelSpace);

% Load toolbox
addpath('Your BMW Toolbox Path')
BMW('AddPath');

MA_All=cell(1,length(ModelSpace)); % Pre-allocation
for i=1:length(ModelSpace) % loop models
    
    % Specify data
    MA.Data='Your Data Path'; % e.g. E:/HAHAHAHA/Data_Hahaha_1st_continuous.mat
    
    % Specify model
    MA.Model.Model=ModelSpace{i};
    MA.Model.Variants={'ContinuousK'}; % Consider a continuous capacity
    
    % Here shows all the model variants
%     MA.Model.Variants={'ContinuousK', 'ResponseNoise', 'Bias', 'Swap', 'Bias Fluctuation', 'Precision Fluctuation',...
%         'Category (Within-Item)', 'Category (Between-Item)'};

    % Configuration
    MA=Configuration_BMW(MA);
    
    % Estimation
    MA=ModelFit_BMW(MA);
    MA_All{i}=MA; % pack the results into a cell variable
    
end

%% Model Comparison
MC=ModelComparison_BMW(MA_All);

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
