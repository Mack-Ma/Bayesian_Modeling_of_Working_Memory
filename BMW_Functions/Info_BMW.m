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
%   'Fit Options'/'Criteria_MAP'
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
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox: https://github.com/BMW-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Info]=Info_BMW(Q, Data)
if nargin==1
    switch Q.Item
        case 'Criteria_MAP'
            Info={'DIC2','WAIC2','LME_GHM'};
        case 'Criteria_MLE'
            Info={'LLH','AIC','BIC'};
        case 'Fit Options'
            Info.RFXBMS=1;
            Info.UniformPrior=1;
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
                        Info.GAOptions.MaxIter=4000;
                        Info.GAOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                    case 'SA'
                        Info.Algorithm='SA';
                        Info.SAOptions.MaxIter=4000;
                        Info.SAOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                    case 'MADS'
                        Info.Algorithm='MADS';
                        Info.MADSOptions.MaxIter=5000;
                        Info.MADSOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                    case 'BADS'
                        Info.Algorithm='BADS';
                        Info.BADSOptions.MaxIter=5000;
                        Info.BADSOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                    case 'fmincon: sqp'
                        Info.Algorithm='fmincon: sqp';
                        Info.fminconOptions.MaxIter=4000;
                        Info.fminconOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                    case 'fmincon: interior-point'
                        Info.Algorithm='fmincon: interior-point';
                        Info.fminconOptions.MaxIter=4000;
                        Info.fminconOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                    case 'fmincon: active-set'
                        Info.Algorithm='fmincon: active-set';
                        Info.fminconOptions.MaxIter=4000;
                        Info.fminconOptions.StepTolerance=1e-6;
                        Info.Display='iter';
                end
            end
    end
else
    if iscell(Data)
        MaxSS=max(Data{1}.SS);
        Nss=length(unique(Data{1}.SS));
    elseif isstruct(Data)
        MaxSS=max(Data.SS);
        Nss=length(unique(Data.SS));
    else
        error('Sorry, the data file is invalid...')
    end
    switch Q.Model
        case 'Item Limit'
            Info.ParamInfo='Item Limit: K (capacity), kappa_r (motor noise)';
            Info.start=[2, 50];
            Info.ub=[MaxSS+1e-6, 700];
            Info.lb=[0, 0.001];
        case 'Standard Mixture'
            Info.ParamInfo='Standard Mixture: K (capacity), kappa(N) (precision at set size N)';
            Info.start=[2, 50*ones(1,Nss)];
            Info.ub=[MaxSS+1e-6, 700*ones(1,Nss)];
            Info.lb=[0, 0.001+zeros(1,Nss)];
        case 'Slots-plus-Averaging'
            Info.ParamInfo='Slots-plus-Averaging: K (capacity), kappa_1 (unit precision), kappa_r (motor noise)';
            Info.start=[2, 50];
            Info.ub=[MaxSS+1e-6, 700];
            Info.lb=[0, 0.001];
        case 'Equal Precision'
            Info.ParamInfo='Equal Precision: kappa1_bar (precision at set size 1), power (precision decay rate), kappa_r (motor noise)';
            Info.start=[50, 1];
            Info.ub=[700, 10];
            Info.lb=[0.001, 0.001];
        case 'Variable Precision'
            if Nss~=1
                Info.ParamInfo='Variable Precision: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate)';
                Info.start=[50, 50, 1];
                Info.ub=[700, 700, 10];
                Info.lb=[0.001, 0.001, 0];
            else
                Info.ParamInfo='Variable Precision: kappa1_bar (mean precision at set size 1), tau';
                Info.start=[50, 50];
                Info.ub=[700, 700];
                Info.lb=[0.001, 0.001];
            end
        case 'Variable Precision with Capacity'
            if Nss~=1
                Info.ParamInfo='Variable Precision with Capacity: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), K (capacity)';
                Info.start=[50, 50, 1, 2];
                Info.ub=[700, 700, 10, MaxSS+1e-6];
                Info.lb=[0.001, 0.001, 0, 0];
            else
                Info.ParamInfo='Variable Precision with Capacity: kappa1_bar (mean precision at set size 1), tau, kappa_r (motor noise), K (capacity)';
                Info.start=[50, 50, 2];
                Info.ub=[700, 700, MaxSS+1e-6];
                Info.lb=[0.001, 0.001, 0];
            end
        case 'Category-Only'
            Info.ParamInfo='Category-Only: kappa_c (categorical precision)';
            Info.start=50;
            Info.ub=700;
            Info.lb=0.001;
        case 'Category-Only (with Capacity)'
            Info.ParamInfo='Category-Only (with Capacity): kappa_c (categorical precision), K (capacity)';
            Info.start=[50 2];
            Info.ub=[700 MaxSS+1e-6];
            Info.lb=[0.001, 0];
        case 'Fixed-Capacity (Full-Display)'
            Info.ParamInfo='Fixed-Capacity: K(capacity), g(guess)';
            Info.start=[2, 0.5];
            Info.ub=[MaxSS+1e-6, 1];
            Info.lb=[0, 0];
        case 'Fixed-Capacity (Single-Probe)'
            Info.ParamInfo='Fixed-Capacity: K(capacity), g(guess)';
            Info.start=[2, 0.5];
            Info.ub=[MaxSS+1e-6, 1];
            Info.lb=[0, 0];
        case 'Fixed-Capacity (Central-Probe)'
            Info.ParamInfo='Fixed-Capacity: K(capacity), g(guess)';
            Info.start=[2, 0.5];
            Info.ub=[MaxSS+1e-6, 1];
            Info.lb=[0, 0];
    end
    if isfield(Q, 'Variants')
        % Continuous recall models
        if any(strcmp(Q.Variants,'Category (Between-Item)')) && any(strcmp(Q.Model,{'Standard Mixture','Slots-plus-Averaging',...
                'Variable Precision','Variable Precision with Capacity'}))
            if any(strcmp(Q.Variants,'VariableCatWeight'))
                Info.ParamInfo=[Info.ParamInfo, ', kappa_c (categorical precision), p_c (categorical memory rate)'];
                Info.start=[Info.start, 50, 0.5*ones(1,Nss)];
                Info.ub=[Info.ub, 700, ones(1,Nss)];
                Info.lb=[Info.lb, 0.001, zeros(1,Nss)];
            else
                Info.ParamInfo=[Info.ParamInfo, ', kappa_c (categorical precision), p_c (categorical memory rate)'];
                Info.start=[Info.start, 50, 0.5];
                Info.ub=[Info.ub, 700, 1];
                Info.lb=[Info.lb, 0.001, 0];
            end
        end 
        if any(strcmp(Q.Variants,'Category (Within-Item)')) && any(strcmp(Q.Model,{'Slots-plus-Averaging','Variable Precision'}))
            if any(strcmp(Q.Variants,'VariableCatWeight'))
                Info.ParamInfo=[Info.ParamInfo, ', kappa_c (categorical precision), p_c (categorical memory rate)'];
                Info.start=[Info.start, 50, 0.5*ones(1,Nss)];
                Info.ub=[Info.ub, 700, ones(1,Nss)];
                Info.lb=[Info.lb, 0.001, zeros(1,Nss)];
            else
                Info.ParamInfo=[Info.ParamInfo, ', kappa_c (categorical precision), eps_c (categorical memory scaling)'];
                Info.start=[Info.start, 50, 0.5];
                Info.ub=[Info.ub, 700, 1];
                Info.lb=[Info.lb, 0.001, -1];
            end
        end
        if any(strcmp(Q.Variants,'ResponseNoise')) && ~any(strcmp(Q.Model,{'Item Limit'}))
            Info.ParamInfo=[Info.ParamInfo, ', kappa_r(response precision)'];
            Info.start=[Info.start, 50];
            Info.ub=[Info.ub, 700];
            Info.lb=[Info.lb, 0.001];
        end
        if any(strcmp(Q.Variants,'Bias'))
            Info.ParamInfo=[Info.ParamInfo, ', mu(bias)'];
            Info.start=[Info.start, 0];
            Info.ub=[Info.ub, 90];
            Info.lb=[Info.lb, -90];
        end
        if any(strcmp(Q.Variants,'BiasF')) && ~any(strcmp(Q.Variants,{'Category (Between-Item)','Category (Within-Item)'}))
            Info.ParamInfo=[Info.ParamInfo, ', muf(bias fluctuation)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 90];
            Info.lb=[Info.lb, -90];
        end
        if any(strcmp(Q.Variants,'PrecF')) && ~any(strcmp(Q.Variants,{'Category (Between-Item)','Category (Within-Item)'})) && ~any(strcmp(Q.Model,'Item Limit'))
            Info.ParamInfo=[Info.ParamInfo, ', kappa_f(precision fluctuation)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 7];
            Info.lb=[Info.lb, -7];
        end
        if any(strcmp(Q.Variants,'Swap'))
            Info.ParamInfo=[Info.ParamInfo, ', b(swap rate)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 1];
            Info.lb=[Info.lb, 0];
        end
        % Change detection models
        if any(strcmp(Q.Variants,'Lapse'))
            Info.ParamInfo=[Info.ParamInfo, ', a(lapse rate)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 1];
            Info.lb=[Info.lb, 0];
        end
        if any(strcmp(Q.Variants,'Ensemble'))
            Info.ParamInfo=[Info.ParamInfo, ', e(rate of ensemble encoding)'];
            Info.start=[Info.start, 0.1];
            Info.ub=[Info.ub, 1];
            Info.lb=[Info.lb, 0];
        end
    end
end