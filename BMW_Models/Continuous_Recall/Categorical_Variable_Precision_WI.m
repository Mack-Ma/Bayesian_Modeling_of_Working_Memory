%% Categorical Slots-plus-Averaging (Within-Variant)
%
% Define the Categorical Slots-plus-Averaging (Within-Item) model
% ------------
% Output=Categorical_Slots_plus_Averaging_WI(param, Data, Input)
%
% ## Theory ##
% This model assumed that the responses were supported by both categorical
% representations and continuous representations. For continuous memory,
% there are a limited number of resouce slots and the fidelity of each representation
% depends on the number of slots being allocated to the corresponding item.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% K, kappa_1, kappa_r, kappa_c, p_c, (bias, s)
%
% - Data
% Data.error (response-sample), Data.SS (set size), Data.error_c, (category-sample)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt, Data.error_nt_c)
%
% - Input
% Input.Variants
%   cell array, options of model variants
% Input.Output
%   string, choose output mode
%       'Prior', only output prior density
%       'LLH', output log likelihood
%       'LP', output log posterior density
% Input.PDF
%   0/1, output pdf or not. default as 0
%   Valid only when Input.Output=='LLH'
%
% ## Output ##
% Output is conditional to Input.Output & Input.PDF
%
% ## Reference ##
% - Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations in visual working memory.
% Nature, 453(7192), 233.
% - van den Berg, R., Shin, H., Chou, W. C., George, R., & Ma, W. J. (2012).
% "Variability in encoding precision accounts for visual short-term memory limitations".
% Proceedings of the National Academy of Sciences, 109(22), 8780-8785.
% - Bays, P. M., Catalao, R. F., & Husain, M. (2009). "The precision of visual working memory
% is set by allocation of a shared resource." Journal of Vision, 9(10), 7-7.
% - Hardman, K. O., Vergauwe, E., & Ricker, T. J. (2017). Categorical working memory representations
% are used in delayed estimation of continuous colors. Journal of Experimental Psychology:
% Human Perception and Performance, 43(1), 30.
%
% ------------
% Programmed by Ma, Tianye
% Mentored by Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%


function Output=Categorical_Variable_Precision_WI(param, Data, Input)

% Specify parameters
kappa1_bar=param(1); % Precision at set size 1
tau=param(2); % Resource allocation variability
Nparam=2;
SS=Data.SS;
SS_range=unique(SS);
if length(SS_range)~=1
    Nparam=Nparam+1;
    power=param(Nparam); % Resource decay rate
else
    power=0;
end
kappa_c=param(Nparam+1); % Precision of categorical memory
eps_c=param(Nparam+2); % Scaling the weight of categorical memory
Nparam=Nparam+2;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'ResponseNoise'))
    Nparam=Nparam+1;
    kappa_r=param(Nparam); % Response precision
end
if ~any(strcmp(Input.Variants,'Bias'))
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if ~any(strcmp(Input.Variants,'Swap'))
    s=0; % No swap
else
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
end

% Configuration
samples=Data.sample;
responses=Data.response;
categories=Data.category;
error_range=Data.error_range;
if length(error_range)==2
    period=error_range(2)-error_range(1);
else
    period=max(error_range)-min(error_range)+(error_range(2)-error_range(1));
end
errors=CircDist_BMW('Diff',responses,samples,period);
errors_c=CircDist_BMW('Diff',responses,categories,period);
response_range=Data.response_range;
sample_range=Data.sample_range;
category_range=unique(Data.category);
if any(strcmp(Input.Variants,'Swap'))
    samples_nt=Data.sample_nt;
    categories_nt=Data.category_nt;
    errors_nt=CircDist_BMW('Diff',repmat(responses,1,size(samples_nt,2)),samples_nt,period);
    errors_nt_c=CircDist_BMW('Diff',repmat(responses,1,size(samples_nt,2)),categories_nt,period);
end
kappa_max=700;
SampleSeed=2000;
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
    % LH function
    p_error=zeros(length(SS_range),length(error_range));
    p_error_c=zeros(length(SS_range),length(error_range));
    for i_N=1:length(SS_range)
        N=SS_range(i_N);
        
        % MC Sampling
        kappa_bar=kappa1_bar*ones(1,SampleSeed)/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        p_error0=zeros(SampleSeed,length(error_range));
        if any(strcmp(Input.Variants,'ResponseNoise'))
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias;
                conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
                conv_kappa_c=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0));
                p_error0(:,i_error)=besseli(0,conv_kappa)./(2*pi*besseli(0,kappa)*besseli(0,kappa_r));
                p_error_c(i_N,i_error)=besseli(0,conv_kappa_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
            end
        else
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias;
                p_error0(:,i_error)=exp(kappa.*cosd(error0))./(2*pi*besseli(0,kappa));
                p_error_c(i_N,i_error)=exp(kappa_c.*cosd(error0))./(2*pi*besseli(0,kappa_c));
            end
        end
        p_error0=p_error0(~isinf(sum(p_error0,2)),:);
        p_error(i_N,:)=mean(p_error0(~isnan(sum(p_error0,2)),:),1); % Find average across samples
        
        % Normalization
        p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
        p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
        
    end
    
    % calculate the normalization constant
    nc=zeros(length(SS_range),length(sample_range));
    period_sample=sample_range(end)+sample_range(2)-2*sample_range(1);
    for i_N=1:length(SS_range)
        for i_s=1:length(sample_range)
            S=sample_range(i_s);
            % find the current category
            [~,C_ind]=min(abs(CircDist_BMW('Diff',repmat(S,size(category_range,1),1),category_range,period_sample)));
            C=category_range(C_ind);
            e_s=CircDist_BMW('Diff',response_range,ones(size(response_range))*S,period_sample); % continuous error
            e_c=CircDist_BMW('Diff',response_range,ones(size(response_range))*C,period_sample); % categorical error
            s_ind=interp1(error_range,1:length(error_range),e_s);
            c_ind=interp1(error_range,1:length(error_range),e_c);
            nc(i_N,i_s)=sum(p_error(i_N,s_ind).*(p_error_c(i_N,c_ind)).^eps_c);
        end
    end
    
    % Calculate LH
    p_T=zeros(1,length(errors));
    p_NT=zeros(1,length(errors));
    for i=1:length(errors)
        p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i))*(p_error_c(SS_range==SS(i),error_range==errors_c(i))).^eps_c/...
            nc(SS_range==SS(i),sample_range==samples(i));
    end
    if any(strcmp(Input.Variants,'Swap'))
        for i=1:length(errors_nt)
            if SS(i)==1
                p_NT(i)=0;
            else
                for j=1:SS(i)-1
                    p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(i,j)).*(p_error_c(SS_range==SS(i),error_range==errors_nt_c(i,j))).^eps_c./nc(SS_range==SS(i),sample_range==samples_nt(i,j));
                end
            end
        end
    end
    p_LH=p_T+p_NT; % Target + non-target LH
    
    % LLH
    if isfield(Input,'PDF') && Input.PDF==1
        LLH.error=p_error; % PDF
    else
        LLH=-sum(log(p_LH)); % Negative LLH
    end
    
    % Posterior
    LP=-log(Prior)+LLH; % likelihood*prior
    
end

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

% Define prior
function p=prior(param, Input)

% Specify parameters
kappa1_bar=param(1); % Unit resource
% Gamma prior for unit resource
p0(1)=gampdf(kappa1_bar,3,15);
tau=param(2); % Resource allocation variability
% prior for tau
% Note that here we simplified the theoretical conjugate prior of
% gamma distribution for convenience
p0(2)=2.^(-0.05*tau)./tau.^(-2)/48040; % normalized
% check power
if length(SS_range)~=1
    Nparam=3;
    power=param(Nparam); % decay rate
    % gamma prior for power
    p0(Nparam)=gampdf(power,1.5,1);
else
    Nparam=2;
end
kappa_c=param(Nparam+1); % categorical memory precision
% Gamma prior for categorical memory precision
p0(Nparam+1)=gampdf(kappa_c,3,5);
eps_c=param(Nparam+2); % categorical weight
% Gamma prior for the weight of categorical memory
p0(Nparam+2)=gampdf(eps_c, 2, 2);
Nparam=Nparam+2;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'ResponseNoise'))
    Nparam=Nparam+1;
    kappa_r=param(Nparam); % Response precision
    % Gamma prior for response precision 
    p0(Nparam)=gampdf(kappa_r,3,5);
end
if any(strcmp(Input.Variants,'Bias'))
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
    % Gaussian prior for bias
    p0(Nparam)=normpdf(bias, 0, 1);
end
if any(strcmp(Input.Variants,'BiasF'))
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
    % Gaussian prior for the fluctuation of bias
    p0(Nparam)=normpdf(biasF, 0, 5);
end
if any(strcmp(Input.Variants,'PrecF'))
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
    % Gaussian prior for the fluctuation of precision
    p0(Nparam)=normpdf(precF, 0, 1);
end
if any(strcmp(Input.Variants,'Swap'))
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
    % Gaussian prior for the swap rate
    p0(Nparam)=normpdf(s, 0.5, 1);
end

% Construct joint distribution
% Consider independent parameters here
% We think it's generally acceptable for prior definition,
% tho it's usually not the actual case
p=1;
for i=1:Nparam
    p=p*p0(i);
end

end
