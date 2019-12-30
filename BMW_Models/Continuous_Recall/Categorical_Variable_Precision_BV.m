%% Between-Variant Categorical Variable Precision
%
% Log likelihood function of the Categorical Variable Precision (Between-Variant) model
% ------------
% LLH=Categorical_Variable_Precision_BV(param, Data, Input)
%
% ## Theory ##
% This model assumed that the responses were supported by both categorical
% representations and continuous representations. For continuous memory,
% there's no item-based capacity limit and the relationship between memory fidelity
% and set size follows a power-law function. In addition, memory resource should have
% trial-by-trial, item-by-item variation.
% Response follows a Von Mises distribution in which the shape parameter
% follows Gamma distribution
%
% ## Input ##
% - param
% kappa_1, tau, power, kappa_r, kappa_c, p_c (bias, biasF, s)
% - Data
% Data.error (response-sample), Data.SS (set size), Data.error_c, (category-sample)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt, Data.error_nt_c)
% - Input
% Input.Variants (options of Variants), Input.PDF
%
% ## Output ##
% Summed log Likelihood
%
% ## Reference ##
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
% Under the guidance of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function LLH=Categorical_Variable_Precision_BV(param, Data, Input)

% Specify parameters
kappa1_bar=param(1);
tau=param(2);
power=param(3);
kappa_r=param(4);
kappa_c=param(5);
p_c=param(6);
Nparam=6;
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
SS=Data.SS;
SS_range=unique(SS);
if Input.Variants.Swap==1
    errors_nt=Data.error_nt;
    errors_nt_c=Data.error_nt_c;
end
kappa_max=700; % Computational limit
SampleSeed=1000; % Monte Carlo seed

% LH function
if continuous==1
    
    p_error0=zeros(SampleSeed,length(errors));
    p_error_c=zeros(1,length(errors));
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
        % Von Mises distribution convoluted by kappa_r
        conv_kappa=sqrt(kappa(:,SS_range==N).^2+kappa_r^2+2*kappa(:,SS_range==N)*kappa_r.*cosd(error0));
        conv_kappa_c=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0));
        p_error0(:,i_error)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r));
        p_error_c(i_error)=besseli0_fast(conv_kappa_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r));
        
        if Input.Variants.Swap==1
            if N==1
                p_error0_NT(:,i_error)=0;
            else
                p_temp_NT=zeros(SampleSeed,1);
                for i_nt=1:N-1
                    error0_nt=errors_nt(i_nt, i_error)+bias; % Errors with bias
                    error0_nt_c=errors_nt(i_nt, i_error);
                    conv_kappa_c_nt=sqrt(kappa_c.^2+kappa_r^2+2*kappa_c*kappa_r.*cosd(error0_nt_c));
                    conv_kappa_nt=sqrt(kappa(:,SS_range==N).^2+kappa_r^2+2*kappa(:,SS_range==N)*kappa_r.*cosd(error0_nt));
                    p_temp_NT=p_temp_NT+(1-p_c)*besseli0_fast(conv_kappa_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r))+...
                        p_c*(besseli0_fast(conv_kappa_c_nt)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r)));
                end
                p_error0_NT(:,i_error)=p_temp_NT/(N-1);
            end
        end
    end
    p_T=(1-s)*((1-p_c)*mean(p_error0,1)+p_c*p_error_c);
    p_NT=s*mean(p_error0_NT,1);
    p_LH=p_T+p_NT;
else
    p_error=zeros(length(SS_range),length(error_range));
    p_error_c=zeros(length(SS_range),length(error_range));
    bias_cur=bias*ones(1,SampleSeed);
    for i_N=1:length(SS_range)
        N=SS_range(i_N);
        
        % MC Sampling
        rng('shuffle');
        kappa_bar=kappa1_bar*ones(1,SampleSeed)/(N).^power;
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
        p_error(i_N,:)=mean(p_error0,1);
        
        % Normalization
        p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
        p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
        
    end
    
    % Calculate LH
    p_T=zeros(1,length(errors));
    p_NT=zeros(1,length(errors));
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
    if Input.Variants.PrecF==1 || Input.Variants.BiasF==1
        if length(unique(Data.sample_range))~=1
            error('PDF mode requires sample to be unique')
        end
    elseif length(unique(Data.SS))~=1
        error('PDF mode requires set size to be unique')
    else
        LLH=p_error; % PDF
    end
else
    LLH=-sum(log(p_LH)); % Negative LLH
    if abs(LLH)==Inf || isnan(LLH)
        LLH=exp(666); % LLH should be a real value
    end
end

end
