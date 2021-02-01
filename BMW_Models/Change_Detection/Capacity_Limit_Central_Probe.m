%% Fixed-Capacity (Central-Probe)
%
% Log likelihood function of the Fixed-Capacity model (Central-Probe)
% ------------
% ## Theory ##
% This model assumed that there's a fixed discrete item-based capacity so
% there should be a proportion of guess response due to the capacity-limit.
% Note that it is for central-probe change detection paradigm only, any kind 
% of misuse will cause trouble. In addition, Rouder et al.(2008) revealed
% the variation of guess rate according to the change probability, so it'll
% be more appropriate to set the change probability as 0.5 consistently.
%
% ## Input ##
% - param
% K, g, (a, e)
% - input
% input.response(report: change/same), input.cond(trial type: change/same)
% input.SS(set size), input.deriv(derivative sets)
%
% ## Output ##
% Summed log likelihood
%
% ## Reference ##
% - Cowan, N. (2001). "The magical number 4 in short-term memory: 
% A reconsideration of mental storage capacity". Behavioral and Brain Sciences, 24(1), 87-114.
% - Rouder, J. N., Morey, R. D., Cowan, N., Zwilling, C. E., Morey, C. C., & Pratte, M. S. (2008). 
% "An assessment of fixed-capacity models of visual working memory". 
% Proceedings of the National Academy of Sciences, 105(16), 5975-5979.
% - Liesefeld, H. R., Liesefeld, A. M., & Müller, H. J. (2019). 
% "Two good reasons to say ‘change!’–ensemble representations 
% as well as item representations impact standard measures of VWM capacity". 
% British Journal of Psychology, 110(2), 328-356.
%
% ------------
% Programmed by Ma, Tianye
% Under the guidance of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 10/21/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%


function Output=Capacity_Limit_Central_Probe(param, Data, Input)

% Specify parameters
if ~isfield(Input,'Variants') % No variants
    Input.Variants={};
end
K=param(1); % Fixed Capacity/Capacity at 1st set size
p_change=Data.p_change;
p_change_range=unique(p_change);
Npchange=length(p_change_range);
g=param(2:1+Npchange); % Guess Rate
Nparam=1+Npchange;
% Specify Variants
SS=Data.SS;
SS_range=unique(SS);
Nss=length(SS_range);
if any(strcmp(Input.Variants,'VariableK')) && Nss>1
    K_others=param(Nparam+1,Nparam+Nss-1); % Variable Capacity
    K_all=[K K_others];
    Nparam=Nparam+Nss-1;
else
    K_all=K*ones(Nss,1);
end
if ~any(strcmp(Input.Variants,'ContinuousK'))
    K_all=floor(K_all); % Note that capacity is a fixed, discrete value here
end
if any(strcmp(Input.Variants,'Lapse'))
    Nparam=Nparam+1;
    a=param(Nparam); % Lapse Rate
else
    a=0;
end
if any(strcmp(Input.Variants,'Ensemble'))
    Nparam=Nparam+1;
    e=param(Nparam); % Rate of ensemble encoding
else
    e=0;
end
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Data, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

% Load data
response=Data.response; % Binary responses
res_range=unique(response);
response=(response-min(res_range))/(max(res_range)-min(res_range)); % response=0/1
Cond=Data.condition;
cond_range=unique(Cond);
Cond=(response-min(cond_range))/(max(cond_range)-min(cond_range)); % condition=0/1

% Likelihood function
p_hit=zeros(Npchange,Nss);
p_fa=zeros(Npchange,Nss);
for pc=1:Npchange
    g_now=g(pc);
    for i_N=1:length(SS_range)
        K_now=K_all(i_N);
        N=SS_range(i_N);
        Pmem=min(min(1,K_now/N)+e,1); % Capacity limit
        % Hit response
        p_hit(pc,i_N)=g_now;
        % False-alarm
        p_fa(pc,i_N)=(1-a)*(1-Pmem)*g_now+a*g_now;
    end
end

% Calculate LH
p_LH=zeros(1,length(response));
for i=1:length(response)
    switch Cond(i)
        case 0
            p_fa_now=p_fa(p_change_range==p_change(i),SS_range==SS(i));
            switch response(i)
                case 0
                    p_LH(i)=1-p_fa_now;
                case 1
                    p_LH(i)=p_fa_now;
            end
        case 1
            p_hit_now=p_hit(p_change_range==p_change(i),SS_range==SS(i));
            switch response(i)
                case 0
                    p_LH(i)=1-p_hit_now;
                case 1
                    p_LH(i)=p_hit_now;
            end
    end
end

% LLH
LLH=-sum(log(p_LH));

% Posterior
LP=-log(Prior)+LLH; % likelihood*prior

% Decide output
if strcmp(Input.Output,'LP')
    Output=LP;
elseif strcmp(Input.Output,'LLH')
    Output=LLH;
elseif strcmp(Input.Output,'Prior')
    Output=Prior;
elseif strcmp(Input.Output,'LPPD')
    Output=log(p_LH);
elseif strcmp(Input.Output,'All')
    Output.LP=LP;
    Output.LLH=LLH;
    Output.Prior=Prior;
    Output.LPPD=log(p_LH);
end

if ~isstruct(Output) && (any(abs(Output))==Inf || any(isnan(Output)))
    Output=realmax('double'); % Output should be a real value
end

end

function p=prior(param, Data, Input)

% Specify parameters
SS=Data.SS;
SS_range=unique(SS);
Nss=length(SS_range);
if any(strcmp(Input.Variants,'VariableK'))
    K=param(1:Nss);
    p0(1:Nss)=wblpdf(K,3.5,3);
else
    K=param(1); % Capacity
    % weibull prior for capacity
    p0(1)=wblpdf(K,3.5,3); % given that K is ofter 3~4
end
p_change=Data.p_change;
p_change_range=unique(p_change);
Npchange=length(p_change_range);
g=param(Nss:Nss+Npchange); % Guess Rate
% Gaussian prior for the swap rate
p0(Nss:Nss+Npchange)=normpdf(g, 0.5, 1);
Nparam=Npchange+1;
% Specify Variants
if any(strcmp(Input.Variants,'Lapse'))
    Nparam=Nparam+1;
    a=param(Nparam); % Lapse Rate
    p0(Nparam)=normpdf(a, 0.5, 1);
end
if any(strcmp(Input.Variants,'Ensemble'))
    Nparam=Nparam+1;
    e=param(Nparam); % Rate of ensemble encoding
    p0(Nparam)=normpdf(e, 0.5, 1);
end

% Construct joint distribution
% Consider independent parameters here
% We think it's generally acceptable for prior definition
p=1;
for i=1:Nparam
    p=p*p0(i);
end

end