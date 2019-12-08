%% Categorical Variable Precision with Capacity (Between-Variant)
%
% Log Likelihood function of the Categorical Variable Precision with Capacity (Between-Variant) model
% ------------
% ## Theory ##
% This model assumed that there's fixed, discrete, item-based capacity limit and
% the relationship between memory fidelity and set size follows a power-law
% function.  In addition, memory resource is supposed to have the property of
% trial-by-trial, item-by-item variation.
% The response probability density distribution follows the weighted addition of
% random guess (uniform) and memory response (Von Mises) in which
% the shape parameter follows Gamma distribution.
%
% ## Input ##
% - param
% kappa_1, tau, power, kappa_r, K, kappa_c, p_c (bias, biasF, s)
% - Data
% Data.error (response-sample), Data.SS (set size), Data.error_c, (category-sample)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt, Data.error_nt_c)
% - Input
% Input.Derivatives (options of derivatives)
%
% ## Output ##
% Summed log Likelihood
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
% Under the guidance of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function LLH=Categorical_Variable_Precision_with_Capacity_BV(param, Data, Input)

% parameters
kappa1_bar=param(1);
tau=param(2);
power=param(3);
kappa_r=param(4);
K=param(5);
kappa_c=param(6);
p_c=param(7);
Nparam=7;
if Input.Derivatives.Bias==0
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if Input.Derivatives.BiasF==1
    warning('biasF is not identified in the categorical models.')
end
if Input.Derivatives.PrecF==1
    warning('precF is not identified in the categorical models.')
end
if Input.Derivatives.Swap==0
    s=0; % No swap
else
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
end

% Configuration
errors=Data.error;
error_range=Data.error_range;
SS=Data.SS;
SS_range=unique(SS);
errors_c=Data.error_c;
if Input.Derivatives.Swap==1
    errors_nt=Data.error_nt;
    errors_nt_c=Data.error_nt_c;
end
kappa_max=700; % Computational limit
SampleSeed=1000; % Monte Carlo seed

% LH function
p_error=zeros(length(SS_range),length(error_range),1);
p_error_c=zeros(length(SS_range),length(error_range),1);
bias_cur=bias*ones(1,SampleSeed);
for i_N=1:length(SS_range)
    N=SS_range(i_N);
    if K<N
        % MC Sampling
        rng('shuffle'); % Generate random seed
        kappa_bar=kappa1_bar*ones(1, SampleSeed)/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        p_error0=zeros(SampleSeed,length(error_range),1);
        p_error0_c=zeros(SampleSeed,length(error_range),1);
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur;
            conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
            conv_kappa_c=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0));
            p_error0(:,i_error,1)=(1-K/N)*1/length(error_range)+...
                (K/N)*besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
            p_error0_c(:,i_error,1)=(1-K/N)*1/length(error_range)+...
                (K/N)*besseli0_fast(conv_kappa_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r));
        end
        p_error(i_N,:,:)=mean(p_error0,1);
        p_error_c(i_N,:,:)=mean(p_error0_c,1);
        % Normalization
        p_error(i_N,:,1)=p_error(i_N,:,1)/sum(p_error(i_N,:,1));
        p_error_c(i_N,:,1)=p_error_c(i_N,:,1)/sum(p_error_c(i_N,:,1));
    else
        % MC Sampling
        rng('shuffle'); % Generate random seed
        kappa_bar=kappa1_bar*ones(1, SampleSeed)/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        p_error0=zeros(SampleSeed,length(error_range),1);
        p_error0_c=zeros(SampleSeed,length(error_range),1);
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur;
            conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
            conv_kappa_c=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0));
            p_error0(:,i_error,:)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
            p_error0_c(:,i_error,:)=besseli0_fast(conv_kappa_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r));
        end
        p_error(i_N,:,1)=mean(p_error0,1);
        p_error_c(i_N,:,:)=mean(p_error0_c,1);
        % Normalization
        p_error(i_N,:,1)=p_error(i_N,:,1)/sum(p_error(i_N,:,1));
        p_error_c(i_N,:,1)=p_error_c(i_N,:,1)/sum(p_error_c(i_N,:,1));
    end
end

% Calculate LH
p_T=zeros(1,length(errors));
p_NT=zeros(1,length(errors_nt));
for i=1:length(errors)
    p_T(i)=(1-s)*((1-p_c)*p_error(SS_range==SS(i),error_range==errors(i),1)+...
        p_c*p_error_c(SS_range==SS(i),error_range==errors_c(i),1));
end
if Input.Derivatives.Swap==1
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

% LLH
LLH=-sum(log(p_LH)); % Negative LLH

if LLH==Inf || isnan(LLH)
    LLH=exp(666); % LLH should be a real value
end

end
