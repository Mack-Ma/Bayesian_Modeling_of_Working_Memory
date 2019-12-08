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


function LLH=Fixed_Capacity_Central_Probe(param, input)

% Specify parameters
K=param(1); % Capacity
g=param(2); % Guess Rate
% Derivatives
if input.deriv.lapse==1
    a=param(3); % Lapse Rate
else
    a=0;
end
if input.deriv.ensemble==1
    if length(param)<4
        e=param(3); % Rate of ensemble encoding
    else
        e=param(4);
    end
else
    e=0;
end

% Load data
response=input.response; % Binary responses
res_range=unique(response);
response=(response-min(res_range))/(max(res_range)-min(res_range)); % response=0/1
Cond=input.cond;
cond_range=unique(Cond);
Cond=(response-min(cond_range))/(max(cond_range)-min(cond_range)); % condition=0/1
SS=input.SS;
SS_range=unique(SS);
K=floor(K); % Note that capacity is a fixed, discrete value here

% Likelihood function
p_response=zeros(length(SS_range),2,2);
for i_N=1:length(SS_range)
    N=SS_range(i_N);
    Pmem=min(1,K/N)+e; % Capacity limit
    % Hit response
    p_response(i_N,1,1)=(1-a)*g+a*g;
    p_response(i_N,1,2)=1-p_response(i_N,1,1);
    % False-alarm
    p_response(i_N,2,1)=(1-a)*(1-Pmem)*g+a*g;
    p_response(i_N,2,2)=1-p_response(i_N,2,1);
end

% Calculate LH
p_LH=zeros(1,length(response));
for i=1:length(response)
    p_LH(i)=p_response(SS_range==SS(i),2-Cond(i),2-response(i));
end

% LLH
LLH=-sum(log(p_LH));

end