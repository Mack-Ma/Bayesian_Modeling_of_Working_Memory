%% Variable Precision with Capacity
%
% Log Likelihood function of the Variable Precision with Capacity model
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
% kappa_1, tau, power, kappa_r, K, (bias, biasF, s)
% - Data
% Data.error (response-sample), Data.SS (set size)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt)
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
% - Pratte, M. S., Park, Y. E., Rademaker, R. L., & Tong, F. (2017). "Accounting for stimulus-specific variation
% in precision reveals a discrete capacity limit in visual working memory."
% Journal of Experimental Psychology: Human Perception and Performance, 43(1), 6.
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

function LLH=Variable_Precision_with_Capacity(param, Data, Input)

% parameters
kappa1_bar=param(1);
tau=param(2);
power=param(3);
kappa_r=param(4);
K=param(5);
Nparam=5;
if Input.Derivatives.Bias==0
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if Input.Derivatives.BiasF==0
    biasF=0; % Set bias as a consistent value
else
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
end
if Input.Derivatives.PrecF==0
    precF=0; % Set precision as consistent within each set size
else
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
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
if Input.Derivatives.BiasF==1 || Input.Derivatives.PrecF==1
    sample_range=Data.sample_range;
    samples=Data.sample;
else
    sample_range=1;
end
if Input.Derivatives.Swap==1
    errors_nt=Data.error_nt;
    if Input.Derivatives.BiasF==1 || Input.Derivatives.PrecF==1
        samples_nt=Data.sample_nt;
    end
end
kappa_max=700; % Computational limit
SampleSeed=1000; % Monte Carlo seed

% LH function
kappa0_bar=exp(log(kappa1_bar)*(cosd(4*sample_range)).^precF);
bias_cur=bias+biasF*cosd(4*sample_range-90);
bias_cur=ones(SampleSeed,1)*bias_cur;
p_error=zeros(length(SS_range),length(error_range),length(sample_range));
for i_N=1:length(SS_range)
    N=SS_range(i_N);
    if K<N % Beyond capacity
        % MC Sampling
        rng('shuffle'); % Generate random seed
        kappa_bar=ones(SampleSeed,1)*kappa0_bar/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        p_error0=zeros(SampleSeed,length(error_range),length(sample_range));
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur;
            conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
            p_error0(:,i_error,:)=(1-K/N)*1/length(error_range)+...
                (K/N)*besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
        end
        p_error(i_N,:,:)=mean(p_error0,1);
        % Normalization
        for i_sample=1:length(sample_range)
            p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
        end
    else
        % MC Sampling
        rng('shuffle'); % Generate random seed
        kappa_bar=ones(SampleSeed,1)*kappa0_bar/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        p_error0=zeros(SampleSeed,length(error_range),length(sample_range));
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur;
            conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
            p_error0(:,i_error,:)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
        end
        p_error(i_N,:,:)=mean(p_error0,1);
        % Normalization
        for i_sample=1:length(sample_range)
            p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
        end
    end
end

% Calculate LH
p_T=zeros(1,length(errors));
p_NT=zeros(1,length(errors));
for i=1:length(errors)
    if Input.Derivatives.BiasF==1 || Input.Derivatives.PrecF==1
        p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), sample_range==samples(i));
    else
        p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), 1);
    end
end
if Input.Derivatives.Swap==1
    for i=1:length(errors_nt)
        if SS(i)==1
            p_NT(i)=0;
        else
            for j=1:SS(i)-1
                if Input.Derivatives.BiasF==1 || Input.Derivatives.PrecF==1
                    p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(j,i), sample_range==samples_nt(j,i));
                else
                    p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(j,i), 1);
                end
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
