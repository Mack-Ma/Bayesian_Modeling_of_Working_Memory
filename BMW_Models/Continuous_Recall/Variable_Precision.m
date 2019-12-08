%% Variable Precision
%
% Log likelihood function of the Variable Precision model
% ------------
% ## Theory ##
% This model assumed that there's no item-based capacity limit and
% the relationship between memory fidelity and set size follows a power-law
% function. In addition, memory resource is supposed to have the property of
% trial-by-trial, item-by-item variation.
% Response follows a Von Mises distribution in which the shape parameter
% follows Gamma distribution
%
% ## Input ##
% - param
% kappa_1, tau, power, kappa_r, (bias, biasF, precF, s)
% - Data
% Data.error (response-sample), Data.SS (set size)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt)
% - Input
% Input.Derivatives (options of derivatives)
% Input.PDF (return probability density function or not)
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

function LLH=Variable_Precision(param, Data, Input)

% Specify parameters
kappa1_bar=param(1);
tau=param(2);
Nparam=2;
SS=Data.SS;
SS_range=unique(SS);
if length(SS_range)~=1
    Nparam=Nparam+1;
    power=param(Nparam);
else
    power=0;
end
Nparam=Nparam+1;
kappa_r=param(Nparam);
if ~isfield(Input,'Derivatives') % No derivatives
    Input.Derivatives.Bias=0;
    Input.Derivatives.Swap=0;
    Input.Derivatives.BiasF=0;
    Input.Derivatives.PrecF=0;
end
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
if ischar(error_range) && strcmp(error_range,'continuous')
    continuous=1;
else
    continuous=0;
end
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
SampleSeed=1000*Nparam; % Monte Carlo seed

% LH function
kappa0_bar=exp(log(kappa1_bar)*(cosd(4*sample_range)).^precF);
bias_cur=bias+biasF*cosd(4*sample_range-90);
bias_cur=ones(SampleSeed,1)*bias_cur;
if continuous==1 % Continuous mode
    p_error0=zeros(SampleSeed,length(errors),length(sample_range));
    kappa=zeros(SampleSeed,length(SS_range));
    for i_N=1:length(SS_range)
        % MC Sampling
        rng('shuffle'); % Update random seed
        kappa_bar=ones(SampleSeed,1)*kappa0_bar/(SS_range(i_N)).^power;
        kappa0=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa(:,i_N)=min(kappa0, kappa_max); % Constricted by the max kappa
    end
    for i_error=1:length(errors)
        N=SS(i_error);
        
        rng('shuffle'); % Update random seed
        kappa_bar=ones(SampleSeed,1)*kappa0_bar/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        error0=errors(i_error)+bias_cur; % Errors with bias
        % Von Mises distribution convoluted by kappa_r
        conv_kappa=sqrt(kappa(:,SS_range==N).^2+kappa_r^2+2*kappa(:,SS_range==N)*kappa_r.*cosd(error0));
        p_error0(:,i_error,:)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa(:,SS_range==N))*besseli0_fast(kappa_r));
    end
    p_error=mean(p_error0,1);
else
    p_error=zeros(length(SS_range),length(error_range), length(sample_range));
    for i_N=1:length(SS_range)
        N=SS_range(i_N);
        
        % MC Sampling
        rng('shuffle'); % Update random seed
        kappa_bar=ones(SampleSeed,1)*kappa0_bar/(N).^power;
        kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa=min(kappa, kappa_max); % Constricted by the max kappa
        
        p_error0=zeros(SampleSeed,length(error_range),length(sample_range));
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur; % Errors with bias
            % Von Mises distribution convoluted by kappa_r
            conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
            p_error0(:,i_error,:)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r));
        end
        p_error(i_N,:,:)=mean(p_error0,1); % Find average across samples
        
        % Normalization
        for i_sample=1:length(sample_range)
            p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
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
end

% LLH
if isfield(Input,'PDF') && Input.PDF==1
    if Input.Derivatives.PrecF==1 || Input.Derivatives.BiasF==1
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
    if LLH==Inf || isnan(LLH)
        LLH=exp(666); % LLH should be a real value
    end
end

end
