%% Categorical Variable Precision with Capacity (Between-Variant)
%
% Define the Categorical Variable Precision with Capacity (Between-Variant) model
% ------------
% Output=Categorical_Variable_Precision_with_Capacity_BV(param, Data, Input)
%
% ## Theory ##
% This model assumed that the responses were supported by both categorical
% representations and continuous representations. For continuous memory,
% there's a fixed, discrete, item-based capacity limit and the relationship between
% memory fidelity and set size follows a power-law function.
% In addition, memory resource should have trial-by-trial, item-by-item variation.
% The response probability density distribution follows the weighted addition of
% random guess (uniform) and memory response (Von Mises) in which
% the shape parameter follows Gamma distribution.
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% kappa1_bar, tau, (power,) kappa_r, K, kappa_c, p_c, (bias, s)
%
% - Data
% Data.error (response-sample), Data.SS (set size), Data.error_c, (category-sample)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt, Data.error_nt_c)
%
% - Input
% Input.Variant
%   options of model variants
%       Input.Variants.Bias, 0/1 to decide whether consider representational
%       shift/response bias
%       Input.Variants.Swap, 0/1 to decide whether use swap variants
%       Input.Variants.BiasF, 0/1 to decide whether consider a cosine-shaped
%       fluctuation of representational shift/response bias
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
% - van den Berg, R., Shin, H., Chou, W. C., George, R., & Ma, W. J. (2012).
% "Variability in encoding precision accounts for visual short-term memory limitations".
% Proceedings of the National Academy of Sciences, 109(22), 8780-8785.
% - van den Berg, R., Awh, E., & Ma, W. J. (2014). Factorial comparison of working memory models.
% Psychological review, 121(1), 124.
% - Bays, P. M., Catalao, R. F., & Husain, M. (2009). "The precision of visual working memory
% is set by allocation of a shared resource." Journal of Vision, 9(10), 7-7.
% - Hardman, K. O., Vergauwe, E., & Ricker, T. J. (2017). Categorical working memory representations
% are used in delayed estimation of continuous colors. Journal of Experimental Psychology:
% Human Perception and Performance, 43(1), 30.
%
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function Output=Categorical_Variable_Precision_with_Capacity_BV(param, Data, Input)

% parameters
kappa1_bar=param(1);
tau=param(2);
SS=Data.SS;
SS_range=unique(SS);
Nparam=2;
if length(SS_range)~=1
    Nparam=Nparam+1;
    power=param(Nparam); % resource decay rate
else
    power=0;
end
kappa_r=param(Nparam+1); % response noise
K=param(Nparam+2)
kappa_c=param(Nparam+3); % categorical memory precision
p_c=param(Nparam+4); % categorical weight
Nparam=Nparam+7;
if Input.Variants.Bias==0
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if Input.Variants.BiasF==1
    warning('biasF is not identified in the categorical models.')
end
if Input.Variants.PrecF==1
    warning('precF is not identified in the categorical models.')
end
if Input.Variants.Swap==0
    s=0; % No swap
else
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
end

% Configuration
errors=Data.error;
error_range=Data.error_range;
if length(error_range)==2
    continuous=1;
else
    continuous=0;
end
errors_c=Data.error_c;
if Input.Variants.Swap==1
    errors_nt=Data.error_nt;
    errors_nt_c=Data.error_nt_c;
end
kappa_max=700; % Computational limit
SampleSeed=1000; % Monte Carlo seed
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior')
    Prior=prior(param, Input, SS_range); % get prior
elseif strcmp(Input.Output,'LLH')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
% LH function
if continuous==1
    
    p_error0=zeros(SampleSeed,length(errors));
    p_error0_NT=zeros(SampleSeed,length(errors));
    kappa=zeros(SampleSeed,length(SS_range));
    
    % MC Sampling
    for i_N=1:length(SS_range)
        kappa_bar=ones(SampleSeed,1)*kappa1_bar/(SS_range(i_N)).^power;
        kappa0=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa(:,i_N)=min(kappa0, kappa_max); % Constricted by the max kappa
    end
    
    for i_error=1:length(errors)
        N=SS(i_error);
        error0=errors(i_error)+bias; % Errors with bias
        error0_c=errors_c(i_error);
        % Von Mises distribution convoluted by kappa_r
        conv_kappa=sqrt(kappa(:,SS_range==N).^2+kappa_r^2+2*kappa(:,SS_range==N)*kappa_r.*cosd(error0));
        conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_c));
        if K<N
            p_error0(:,i_error)=(1-K/N)*1/(error_range(2)-error_range(1))+...
                (K/N)*((1-p_c)*besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r))+...
                p_c*besseli0_fast(conv_1_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r)));
        else
            p_error0(:,i_error)=(1-p_c)*besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r))+...
                p_c*besseli0_fast(conv_1_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r));
        end
        
        if Input.Variants.Swap==1
            if N==1
                p_error0_NT(:,i_error)=0;
            else
                p_temp_NT=zeros(SampleSeed,1);
                if K<N
                    for i_nt=1:N-1
                        error0_nt=errors_nt(i_nt, i_error)+bias; % Errors with bias
                        conv_1_c_nt=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_nt_c));
                        conv_kappa_nt=sqrt(kappa(:,SS_range==N).^2+kappa_r^2+2*kappa(:,SS_range==N)*kappa_r.*cosd(error0_nt));
                        p_temp_NT=p_temp_NT+(1-K/N)*1/(error_range(2)-error_range(1))+...
                            (K/N)*((1-p_c)*besseli0_fast(conv_kappa_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r))+...
                            p_c*besseli0_fast(conv_1_c_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r)));
                    end
                else
                    for i_nt=1:N-1
                        error0_nt=errors_nt(i_nt, i_error)+bias; % Errors with bias
                        conv_1_c_nt=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_nt_c));
                        conv_kappa_nt=sqrt(kappa(:,SS_range==N,i_error).^2+kappa_r^2+2*kappa(:,SS_range==N,i_error)*kappa_r.*cosd(error0_nt));
                        p_temp_NT=p_temp_NT+(1-p_c)*besseli0_fast(conv_kappa_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N,i_error))*besseli0_fast(kappa_r))+...
                            p_c*besseli0_fast(conv_1_c_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r));
                    end
                end
                p_error0_NT(:,i_error)=p_temp_NT/(N-1);
            end
        end
    end
    p_error0=p_error0(~isinf(sum(p_error0,2)),:);
    p_T=(1-s)*mean(p_error0(~isnan(sum(p_error0,2)),:),1);
    p_NT=s*mean(p_error0_NT,1);
    p_LH=p_T+p_NT;
    
else
    p_error=zeros(length(SS_range),length(error_range));
    p_error_c=zeros(length(SS_range),length(error_range));
    bias_cur=bias*ones(1,SampleSeed);
    for i_N=1:length(SS_range)
        N=SS_range(i_N);
        if K<N
            % MC Sampling
            kappa_bar=kappa1_bar*ones(1, SampleSeed)/(N).^power;
            kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
            kappa=min(kappa, kappa_max); % Constricted by the max kappa
            
            p_error0=zeros(SampleSeed,length(error_range));
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias_cur;
                conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
                conv_kappa_c=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0));
                p_error0(:,i_error)=(1-K/N)*1/length(error_range)+...
                    (K/N)*besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
                p_error_c(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                    (K/N)*besseli0_fast(conv_kappa_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r));
            end
            p_error0=p_error0(~isinf(sum(p_error0,2)),:);
            p_error(i_N,:)=mean(p_error0(~isnan(sum(p_error0,2)),:),1); % Find average across samples
            
            % Normalization
            p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
            p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
        else
            % MC Sampling
            rng('shuffle'); % Generate random seed
            kappa_bar=kappa1_bar*ones(1, SampleSeed)/(N).^power;
            kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
            kappa=min(kappa, kappa_max); % Constricted by the max kappa
            
            p_error0=zeros(SampleSeed,length(error_range));
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias_cur;
                conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
                conv_kappa_c=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0));
                p_error0(:,i_error)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
                p_error_c(i_N,i_error)=besseli0_fast(conv_kappa_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r));
            end
            p_error0=p_error0(~isinf(sum(p_error0,2)),:);
            p_error(i_N,:)=mean(p_error0(~isnan(sum(p_error0,2)),:),1); % Find average across samples
            % Normalization
            p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
            p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
        end
    end
    
    % Calculate LH
    p_T=zeros(1,length(errors));
    p_NT=zeros(1,length(errors_nt));
    for i=1:length(errors)
        p_T(i)=(1-s)*((1-p_c)*p_error(SS_range==SS(i),error_range==errors(i),1)+...
            p_c*p_error_c(SS_range==SS(i),error_range==errors_c(i),1));
    end
    if Input.Variants.Swap==1
        for i=1:length(errors_nt)
            if SS(i)==1
                p_NT(i)=0;
            else
                for j=1:SS(i)-1
                    p_NT(i)=p_NT(i)+s/(SS(i)-1)*((1-p_c)*p_error(SS_range==SS(i),error_range==errors_nt(j,i), 1)+...
                        p_c*p_error_c(SS_range==SS(i),error_range==errors_nt_c(j,i), 1));
                end
            end
        end
    end
    p_LH=p_T+p_NT; % Target + non-target LH
end

% LLH
if isfield(Input,'PDF') && Input.PDF==1
    LLH=p_error; % PDF
else
    LLH=-sum(log(p_LH)); % Negative LLH
    if abs(LLH)==Inf || isnan(LLH)
        LLH=exp(666); % LLH should be a real value
    end
end

    % Posterior
    LP=log(Prior)+LLH; % likelihood*prior
    
end

% Decide output
if strcmp(Input.Output,'LP')
    Output=LP;
elseif strcmp(Input.Output,'LLH')
    Output=LLH;
elseif strcmp(Input.Output,'Prior')
    Output=Prior;
end

end

% Define prior
function p=prior(param, Input, SS_range)

% Specify parameters
kappa1_bar=param(1); % Unit resource
% Gamma prior for unit resource
p0(1)=gampdf(kappa1_bar,3,15);
tau=param(2); % Resource allocation variability
% prior for tau
% Note that here we simplified the theoretical conjugate prior of
% gamma distribution for convenience
p0(2)=2.^(-0.05*tau)./tau.^(-2);
% check power
if length(SS_range)~=1
    Nparam=3;
    power=param(Nparam); % decay rate
    % exponential prior for power (improper)
    p0(Nparam)=exp(-power);
else
    Nparam=2;
end
kappa_r=param(Nparam+1); % Response variability
% Cauchy prior for response noise
% copy-paste from MemToolbox here
p0(Nparam+1)=2./(pi+kappa_r.^2);
K=param(Nparam+2); % Capacity
% weibull prior for capacity
p0(Nparam+2)=wblpdf(K,3.5,3); % given that K is ofter 3~4
kappa_c=param(Nparam+3); % categorical memory precision
% Gamma prior for categorical memory precision
p0(Nparam+3)=gampdf(kappa_c,3,5);
p_c=param(Nparam+4); % categorical weight
% Gaussian prior for the rate of categorical memory
p0(Nparam+4)=normpdf(p_c, 0.5, 1);
Nparam=Nparam+4;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants.Bias=0;
    Input.Variants.Swap=0;
end
if Input.Variants.Bias==1
    bias=param(Nparam); % Mean bias
    % Gaussian prior for bias
    p0(Nparam)=normpdf(bias, 0, 1);
end
if Input.Variants.Swap==1
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
