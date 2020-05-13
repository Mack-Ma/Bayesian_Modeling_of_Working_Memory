%% Configuration_BMW
%
% Configure model estimation & assessment
% -----------------------
% [MA_BMW]=Configuration_BMW(Input)
%
% ## Input ##
% For detailed information, plz refer to the manual (type BMW('Manual') in the command line)
% For the complete example, plz refer to "SampleScript1" in the toolbox
% Input
%   ~.Model
%       struct, define the model of interest
%       ~.Model
%           string, define the name of the model to be estimated
%           This branch is NECESSARY, which means it needs to be defined
%           manually.
%           e.g. MA.Model.Model='Categorical Variable Precision';
%       ~.Output
%           string, 'LP'/'LLH', define the output mode of the models
%           This will determine the general method of model fitting (MAP/MLE) 
%           This branch is optional, the default value is 'LP'.
%           e.g. MA.Model.Output='LLH';
%       ~.Variants
%           struct, set the model variants
%           assign the variant(s) of interest as one
%           This is optional, the default set uses no variants.
%           e.g. MA.Model.Variants.Bias=1;
%   ~.Data
%       string/cell
%       This is NECESSARY.
%       If it is a string, then it should define the directory of data.
%       e.g. MA.Model.Data='E:/matlab/Data_BMW.mat';
%       If it's a cell variable. then it should be the data variable
%       e.g. MA.Model.Data=Data_BMW;
%   ~.FitOptions
%       struct, change the preferences for model fitting
%       Resetting any subfield is optional.
%       ~.Algorithm
%           string, choose the optimization algorithm
%           e.g. MA.FitOptions.Algorithm='GA';
%       ~.UniformPrior
%           boolean variable, choose whether to use the non-informative prior
%           Note that this is only valid for MAP.
%           e.g. MA.FitOptions.UniformPrior=1; % use non-informative prior
%       For other options, plz refer to the manual.
%   ~.Criteria
%       cell, define the criteria of interest
%       This is optional, the default set only uses DIC2, WAIC2 & LME_GHM
%       for MAP and LLH, AIC & BIC for MLE
%       e.g. MA.Criteria={'WAIC1','WAIC2','LME_BS'};
%   ~.Constraints
%       vector, define the initial guess and the constraints of the parameters
%       This is optional.
%       e.g. MA.Constraints.start=[3,50,50];
%
% ## Output ##
% MA_BMW (feed this file to ModelFit_BMW)
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% Sun Yat-Sen University
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [MA_BMW]=Configuration_BMW(Input)

fprintf('\n---------------------------------------- \n')
fprintf('\nNow start configuration... \n')

if ~isfield(Input,'Model')
    error('Sorry, the model of interest should be defined...')
else
    MA_BMW.Model=Input.Model;
end
if ischar(Input.Data)
    RawData=load(Input.Data);
    DataName=fieldnames(RawData);
    if length(DataName)>1
        error('Sorry, the data file is not valid...')
    end
    eval(['MA_BMW.Data=RawData.',DataName{1},';'])
elseif iscell(Input.Data)
    MA_BMW.Data=Input.Data;
else
    error('Sorry, the data file is not valid...')
end
if ~iscell(MA_BMW.Data)
    MA_BMW.Data={MA_BMW.Data};
end

% Default set
if ~isfield(Input,'FitOptions')
    Q.Item='Fit Options'; Q.Algorithm='Default';
    Input.FitOptions=Info_BMW(Q); 
end
if ~isfield(Input.FitOptions,'Method')
    Input.FitOptions.Method='MAP'; % set MAP as default
end
if ~isfield(Input,'Criteria')
    if strcmp(Input.FitOptions.Method,'MAP')
        Q.Item='Criteria_MAP';
        Input.Criteria=Info_BMW(Q);
    elseif strcmp(Input.FitOptions.Method,'MLE')
        Q.Item='Criteria_MLE';
        Input.Criteria=Info_BMW(Q);
    end
end
if ~isfield(Input.Model,'Variants')
    MA_BMW.Model.Variants={}; % No variants
end
if ~isfield(Input,'Constraints')
    QModel.Model=MA_BMW.Model.Model;
    QModel.Variants=MA_BMW.Model.Variants;
    Input.Constraints=Info_BMW(QModel, MA_BMW.Data);
end

% Output
MA_BMW.FitOptions=ConfigFitOptions_BMW(Input.FitOptions);
MA_BMW.Criteria=Input.Criteria;
MA_BMW.Constraints=Input.Constraints;
fprintf('\nDone: Configuration \n')

end