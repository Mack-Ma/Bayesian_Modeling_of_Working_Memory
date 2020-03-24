%% Info_BMW
%
% Return the default set and basic information
% e.g. the default start values of the parameters, the default set of fit
% options etc.
% ------------
% Info=Info_BMW(Q, Data)
%
% ## Input ##
% - Q
% Structual array designate the information of interest
% Q.Item
%   string, designate the information of interest
%   Valid only when Data is not defined.
%   'Fit Options'/'Criteria'
% Q.Algorithm
%   string, if we use 'Fit Options', then the field "Algorithm" will be 
%   needed as well. 
%   e.g. Q.Algorithm='DE-MCMC';
% Q.Model & Q.Variants
%   Designate the model of interest
% e.g. Q.Model='Item Limit', Q.Variants.Swap=1;
%
% - Data: Data file for parameter estimation.
% For the structure of Data, please refer to ModelDefinition_BMW.m
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% Sun Yat-Sen University
% 10/5/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Info]=Info_BMW(Q, Data)
if nargin==1
    switch Q.Item
        case 'Criteria'
            Info.LLH=0;
            Info.AIC=0;
            Info.AICc=0;
            Info.BIC=0;
            Info.DIC=1;
            Info.DICs=0;
            Info.WAIC1=0;
            Info.WAIC2=1;
            Info.LME=1;
        case 'Fit Options'
            if ~isfield(Q,'Algorithm')
                error('Sorry, the field "Algorithm" is needed...')
            else
                switch Q.Algorithm
                    case 'Default'
                        Info.Algorithm='DE-MCMC';
                        Info.Display='iter';
                    case 'DE-MCMC'
                        Info.Algorithm='DE-MCMC';
                        Info.Display='iter';
                    case 'MH-MCMC'
                        Info.Algorithm='MH-MCMC';
                        Info.Display='iter';
                    case 'GA'
                        Info.Algorithm='GA';
                        Info.MaxIter=4000;
                        Info.Display='iter';
                    case 'SA'
                        Info.Algorithm='SA';
                        Info.MaxIter=4000;
                        Info.Display='iter';
                    case 'MADS'
                        Info.Algorithm='MADS';
                        Info.MaxIter=5000;
                        Info.Display='iter';
                    case 'BADS'
                        Info.Algorithm='BADS';
                        Info.MaxIter=5000;
                        Info.Display='iter';
                    case 'fmincon: sqp'
                        Info.Algorithm='fmincon: sqp';
                        Info.MaxIter=4000;
                        Info.Display='iter';
                    case 'fmincon: interior-point'
                        Info.Algorithm='fmincon: interior-point';
                        Info.MaxIter=4000;
                        Info.Display='iter';
                     case 'fmincon: active-set'
                        Info.Algorithm='fmincon: active-set';
                        Info.MaxIter=4000;
                        Info.Display='iter';
                end
            end
    end
else
    if length(Data)>1
        Nss=length(unique(Data{1}.SS));
    elseif length(Data)==1
        Nss=length(unique(Data.SS));
    end
    switch Q.Model
        case 'Item Limit'
            Info.ParamInfo='Item Limit: K (capacity), kappa_r (motor noise)';
            Info.start=[2, 50];
            Info.ub=[20, 700];
            Info.lb=[0, 0.001];
        case 'Standard Mixture'
            Info.ParamInfo='Standard Mixture: K (capacity), kappa(N) (precision at set size N)';
            Info.start=[2, 50*ones(1,Nss)];
            Info.ub=[20, 700*ones(1,Nss)];
            Info.lb=[0, 0.001+zeros(1,Nss)];
        case 'Slots-plus-Averaging'
            Info.ParamInfo='Slots-plus-Averaging: K (capacity), kappa_1 (unit precision), kappa_r (motor noise)';
            Info.start=[2, 50, 50];
            Info.ub=[20, 700, 700];
            Info.lb=[0, 0.001, 0.001];
        case 'Equal Precision'
            Info.ParamInfo='Equal Precision: kappa1_bar (precision at set size 1), power (precision decay rate), kappa_r (motor noise)';
            Info.start=[50, 1, 50];
            Info.ub=[700, 7, 700];
            Info.lb=[0.001, 0.001, 0.001];
        case 'Variable Precision'
            Info.ParamInfo='Variable Precision: kappa1_bar (mean precision at set size 1), tau, power , kappa_r (motor noise)';
            Info.start=[50, 1, 1, 50];
            Info.ub=[700, 500, 10, 700];
            Info.lb=[0.001, 0.001, 0, 0.001];
        case 'Variable Precision with Capacity'
            Info.ParamInfo='Variable Precision with Capacity: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), K (capacity)';
            Info.start=[50, 1, 1, 50, 2];
            Info.ub=[700, 500, 10, 700, 20];
            Info.lb=[0.001, 0.001, 0, 0.001, 0];
        case 'Categorical Slots-plus-Averaging (Between-Variant)'
            Info.ParamInfo='Categorical Slots-plus-Averaging: K (capacity), kappa_1 (unit precision), kappa_r (motor noise), kappa_c(categorical precision), p_c(categorical memory rate)';
            Info.start=[2, 50, 50, 50, .3];
            Info.ub=[20, 700, 700, 700, 1];
            Info.lb=[0, 0.001, 0.001, 0.001, 0];
        case 'Categorical Slots-plus-Averaging (Within-Variant)'
            Info.ParamInfo='Categorical Slots-plus-Averaging: K (capacity), kappa_1 (unit precision), kappa_r (motor noise), kappa_c(categorical precision), p_c(categorical memory rate)';
            Info.start=[2, 50, 50, 50];
            Info.ub=[20, 700, 700, 700];
            Info.lb=[0, 0.001, 0.001, 0.001];
        case 'Categorical Variable Precision (Between-Variant)'
            Info.ParamInfo='Categorical Variable Precision: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), kappa_c (categorical precision), p_c(categorical memory rate)';
            Info.start=[50, 1, 1, 50, 50, .3];
            Info.ub=[700, 500, 10, 700, 700, 1];
            Info.lb=[0.001, 0.001, 0, 0.001, 0.001, 0];
        case 'Categorical Variable Precision with Capacity (Between-Variant)'
            Info.ParamInfo='Categorical Variable Precision with Capacity: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), K (capacity), kappa_c (categorical precision), p_c(categorical memory rate)';
            Info.start=[50, 1, 1, 50, 2, 50, .3];
            Info.ub=[700, 500, 10, 700, 20, 700, 1];
            Info.lb=[0.001, 0.001, 0, 0.001, 0, 0.001, 0];
        case 'Fixed-Capacity (Full-Display)'
            Info.ParamInfo='Fixed-Capacity: K(capacity), g(guess)';
            Info.start=[2, 0.5];
            Info.ub=[20, 1];
            Info.lb=[0, 0];
        case 'Fixed-Capacity (Single-Probe)'
            Info.ParamInfo='Fixed-Capacity: K(capacity), g(guess)';
            Info.start=[2, 0.5];
            Info.ub=[20, 1];
            Info.lb=[0, 0];
        case 'Fixed-Capacity (Central-Probe)'
            Info.ParamInfo='Fixed-Capacity: K(capacity), g(guess)';
            Info.start=[2, 0.5];
            Info.ub=[20, 1];
            Info.lb=[0, 0];
    end
    if isfield(Q, 'Variants')
        % Continuous recall models
        if isfield(Q.Variants, 'Bias') && Q.Variants.Bias==1
            Info.ParamInfo=[Info.ParamInfo, ', mu(bias)'];
            Info.start=[Info.start, 0];
            Info.ub=[Info.ub, 90];
            Info.lb=[Info.lb, -90];
        end
        if isfield(Q.Variants, 'BiasF') && Q.Variants.BiasF==1
            Info.ParamInfo=[Info.ParamInfo, ', muf(bias fluctuation)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 90];
            Info.lb=[Info.lb, -90];
        end
        if isfield(Q.Variants, 'PrecF') && Q.Variants.PrecF==1
            Info.ParamInfo=[Info.ParamInfo, ', kappa_f(precision fluctuation)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 7];
            Info.lb=[Info.lb, -7];
        end
        if isfield(Q.Variants, 'Swap') && Q.Variants.Swap==1
            Info.ParamInfo=[Info.ParamInfo, ', b(swap rate)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 1];
            Info.lb=[Info.lb, 0];
        end
        % Change detection models
        if isfield(Q.Variants, 'Lapse') && Q.Variants.Lapse==1
            Info.ParamInfo=[Info.ParamInfo, ', a(lapse rate)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 1];
            Info.lb=[Info.lb, 0];
        end
        if isfield(Q.Variants, 'Ensemble') && Q.Variants.Ensemble==1
            Info.ParamInfo=[Info.ParamInfo, ', e(rate of ensemble encoding)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 1];
            Info.lb=[Info.lb, 0];
        end
    end
end