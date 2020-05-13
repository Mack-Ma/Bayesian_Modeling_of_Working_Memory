%% Bayesian Modeling of Working Memory (BMW) Toolbox
%
% Fit and compare models of working memory
% -----------------------
% ## Experimental Paradigm ##
% Continuous Recall/ Change Detection/ Custom
%
% ## Model space ##
% ### Continuous Recall ###
% Standard Mixture(Mix)/Item Limit(IL)/Slots-plus-Averaging(SA)/
% Variable Precision(VP)/VP with capacity(VPcap)/
% Equal Precision(EP)/Categorical SA(cSA)/Categorical VP(cVP)/
% Categorical Mixture(cMix)/Categorical VPcap(cVPcap)
% Model Variants: Swap Rate/Bias/Precision Fluctuation/Bias Fluctuation
% ### Change Detection ###
% Fixed-Capacity(Central-Probe/Full-Display/Single-Probe)/Signal Detection
% Model Variants: Lapse Rate/Ensemble Encoding
%
% ## Criteria for model assessment ##
% Logarithmic Likelihood(LLH)/Akaike Information Criterion(AIC)/
% Modified Akaike Information Criterion(AICc)/Bayesian Information Criterion(BIC)/
% Deviance Information Criterion(DIC1/DIC2/DIC*)/
% Watanabe-Akaike Information Criterion(WAIC1/WAIC2)/
% Log Model Evidence Based on the Generalized Harmonic Mean Estimator(LME_GHM)/
% Log Model Evidence Based on the Bridge Sampling Estimator(LME_BS)/
% 2nd-Level Model Frequency/2nd-Level Exceedance Probability
%
% ## Fitting method ##
% Maximum Likelihood(MLE)/Maximum A Posteriori(MAP)
%
% ## Optimization algorithm ##
% fmincon/bads/mads/genetic algorithm/simulated annealing/
% MH-MCMC/DE-MCMC(default)
% ------------
% Programmed by Ma, Tianye
% Mentored by Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% Sun Yat-Sen University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function BMW(varargin)
if nargin==0
    clc % Clean the command window
    fprintf('Howdy!\n')
    fprintf('Bayesian Modeling of Working Memory (BMW) Toolbox\n')
    fprintf('URL: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory\n')
    fprintf('Start GUI...\n')
    AddAllPath_BMW;
    BMW_Main;
    fprintf('Type ''''BMW(''CheckMethods'')'''' to check your available optimization algorithms.\n')
else
    switch varargin{1}
        case 'ClearPath'
            RemoveAllPath_BMW;
        case 'Manual'
            open('BMW_Manual.pdf')
        case 'AddPath'
            AddAllPath_BMW;
        case 'CheckMethods'
            OptList={'DE-MCMC (default)','MH-MCMC','fmincon: sqp','fmincon: active-set','fmincon: interior-point','GA','SA','MADS','BADS'};
            functionList={'BMW_MCMC','BMW_MCMC','fmincon','fmincon','fmincon','ga','simulannealbnd','patternsearch','bads'};
            OKList_ind=[];
            for i=1:length(OptList)
                if exist([functionList{i},'.m'],'file')==2
                    OKList_ind=[OKList_ind,i];                    
                end
            end
            OKList=OptList(OKList_ind);
            fprintf('\nAvailable optimization algorithms:\n\n')
            for i=1:length(OKList)
                fprintf([OKList{i},'\n'])
            end
            fprintf('\r')
    end
end
end

function AddAllPath_BMW
% search BMW.m
path_BMW=which('BMW');
path_BMW=path_BMW(1:end-5);
% add all subpaths
path_list={'BMW_GUI','BMW_Functions','BMW_Models'};
for i=1:3
    path_cur=[path_BMW, path_list{i}];
    addpath(path_cur)
    path_raw=genpath(path_cur);
    seg_ind=[0,find(path_cur==';')];
    if ~isempty(seg_ind)
        for j=1:length(seg_ind)-1
            path_rm=path_raw(seg_ind(j-1)+1:seg_ind(j)-1);
            addpath(path_rm);
        end
    end
end
end

function RemoveAllPath_BMW
% search BMW.m
path_BMW=which('BMW');
path_BMW=path_BMW(1:end-5);
% remove all subpaths
path_list={'BMW_GUI','BMW_Functions','BMW_Models'};
for i=1:3
    path_cur=[path_BMW, path_list{i}];
    rmpath(path_cur)
    path_raw=genpath(path_cur);
    seg_ind=[0,find(path_cur==';')];
    if ~isempty(seg_ind)
        for j=1:length(seg_ind)-1
            path_rm=path_raw(seg_ind(j-1)+1:seg_ind(j)-1);
            rmpath(path_rm);
        end
    end
end
end
