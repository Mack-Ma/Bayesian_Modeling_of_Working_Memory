%% Info_BMW
%
% Return the default set and basic information of the designated model
% ------------
% Info=Info_BMW(Q, Data)
% ## Input ##
% - Q: Designate the issue of interest
% e.g. 'Criteria'/'Fit Options'
% When Q is assigned as a structural array, it should represent the model
% of interest. 
% e.g. Q.Model='Item Limit', Q.Derivatives.Swap=1;
% - Data: Data file recruited to estimated the model.
% For the structure of Data, please refer to ModelDefinition_BMW.m 
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/5/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Info]=Info_BMW(Q, Data)
if nargin==1
    switch Q
        case 'Criteria'
            Info.LLH=1;
            Info.AIC=1;
            Info.AICc=1;
            Info.BIC=1;
            Info.DIC=0;
            Info.WAIC=0;
            Info.LME=0;
        case 'Fit Options'
            Info.Algorithm='GA';
            Info.MaxIter=3000;
            Info.Display='iter';
    end
else
    Nss=length(unique(Data{1}.SS));
    switch Q.Model
        case 'Item Limit'
            Info.ParamInfo='Item Limit: K (capacity), kappa_r (motor noise)';
            Info.start=[2, 50];
            Info.ub=[20, 700];
            Info.lb=[0, 0];
        case 'Standard Mixture'
            Info.ParamInfo='Standard Mixture: K (capacity), kappa(N) (precision at set size N)';
            Info.start=[2, 50*ones(1,Nss)];
            Info.ub=[20, 700*ones(1,Nss)];
            Info.lb=[0, zeros(1,Nss)];
        case 'Slots-plus-Averaging'
            Info.ParamInfo='Slots-plus-Averaging: K (capacity), kappa_1 (unit precision), kappa_r (motor noise)';
            Info.start=[2, 50, 50];
            Info.ub=[20, 700, 700];
            Info.lb=[0, 0, 0];
        case 'Equal Precision'
            Info.ParamInfo='Equal Precision: kappa1_bar (precision at set size 1), power (precision decay rate), kappa_r (motor noise)';
            Info.start=[50, 1, 50];
            Info.ub=[700, 7, 700];
            Info.lb=[0, 0, 0];
        case 'Variable Precision'
            Info.ParamInfo='Variable Precision: kappa1_bar (mean precision at set size 1), tau, power , kappa_r (motor noise)';
            Info.start=[50, 1, 1, 50];
            Info.ub=[700, 500, 10, 700];
            Info.lb=[0, 0, 0, 0];
        case 'Variable Precision with Capacity'
            Info.ParamInfo='Variable Precision with Capacity: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), K (capacity)';
            Info.start=[50, 1, 1, 50, 2];
            Info.ub=[700, 500, 10, 700, 20];
            Info.lb=[0, 0, 0, 0, 0];
        case 'Categorical Slots-plus-Averaging'
            Info.ParamInfo='Categorical Slots-plus-Averaging: K (capacity), kappa_1 (unit precision), kappa_r (motor noise), kappa_c(categorical precision), p_c(categorical memory rate)';
            Info.start=[2, 50, 50, 50, .3];
            Info.ub=[20, 700, 700, 700, 1];
            Info.lb=[0, 0, 0, 0, 0];
        case 'Categorical Variable Precision'
            Info.ParamInfo='Categorical Variable Precision: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), kappa_c (categorical precision), p_c(categorical memory rate)';
            Info.start=[50, 1, 1, 50, 50, .3];
            Info.ub=[700, 500, 10, 700, 700, 1];
            Info.lb=[0, 0, 0, 0, 0, 0];
        case 'Categorical Variable Precision with Capacity'
            Info.ParamInfo='Categorical Variable Precision with Capacity: kappa1_bar (mean precision at set size 1), tau, power (precision decay rate), kappa_r (motor noise), K (capacity), kappa_c (categorical precision), p_c(categorical memory rate)';
            Info.start=[50, 1, 1, 50, 2, 50, .3];
            Info.ub=[700, 500, 10, 700, 20, 700, 1];
            Info.lb=[0, 0, 0, 0, 0, 0, 0];
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
    if isfield(Q.Derivatives, 'Bias')
        switch Q.Derivatives.Bias
            case 1
                Info.ParamInfo=[Info.ParamInfo, ', mu(bias)'];
                Info.start=[Info.start, 0];
                Info.ub=[Info.ub, 90];
                Info.lb=[Info.lb, -90];
        end
        switch Q.Derivatives.BiasF
            case 1
                Info.ParamInfo=[Info.ParamInfo, ', muf(bias fluctuation)'];
                Info.start=[Info.start, 0];
                Info.ub=[Info.ub, 90];
                Info.lb=[Info.lb, -90];
        end
        switch Q.Derivatives.PrecF
            case 1
                Info.ParamInfo=[Info.ParamInfo, ', kappa_f(precision fluctuation)'];
                Info.start=[Info.start, 0];
                Info.ub=[Info.ub, 6];
                Info.lb=[Info.lb, -6];
        end
        switch Q.Derivatives.Swap
            case 1
                Info.ParamInfo=[Info.ParamInfo, ', b(swap rate)'];
                Info.start=[Info.start, 0.1];
                Info.ub=[Info.ub, 1];
                Info.lb=[Info.lb, 0];
        end
    end
    if isfield(Q.Derivatives, 'Lapse') || isfield(Q.Derivatives, 'Ensemble')
        switch Q.Derivatives.Lapse
            case 1
                Info.ParamInfo=[Info.ParamInfo, ', a(lapse rate)'];
                Info.start=[Info.start, 0.1];
                Info.ub=[Info.ub, 1];
                Info.lb=[Info.lb, 0];
        end
        switch Q.Derivatives.Ensemble
            case 1
                Info.ParamInfo=[Info.ParamInfo, ', e(rate of ensemble encoding)'];
                Info.start=[Info.start, 0.1];
                Info.ub=[Info.ub, 1];
                Info.lb=[Info.lb, 0];
        end
    end
end