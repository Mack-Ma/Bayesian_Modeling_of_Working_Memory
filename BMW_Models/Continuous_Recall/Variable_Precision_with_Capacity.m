%% Variable Precision with Capacity
%
% Define the Variable Precision with Capacity model
% ------------
% Output=Variable_Precision_with_Capacity(param, Data, Input)
%
% ## Theory ##
% This model assumed that there's a fixed, discrete, item-based capacity limit and
% the relationship between memory fidelity and set size follows a power-law
% function. In addition, memory resource should have trial-by-trial, item-by-item variation.
% The response probability density distribution follows the weighted addition of
% random guess (uniform) and memory response (Von Mises) in which
% the shape parameter follows Gamma distribution.
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% kappa1_bar, tau, (power,) kappa_r, K, (bias, biasF, precF, s)
%
% - Data
% Data.error (response-sample), Data.SS (set size)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt)
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

function Output=Variable_Precision_with_Capacity(param, Data, Input)

kappa1_bar=param(1); % Precision at set size 1
tau=param(2); % Resource allocation variability
if tau==0, tau=1e-6; end;
Nparam=2;
SS=Data.SS;
SS_range=unique(SS);
if length(SS_range)~=1
    Nparam=Nparam+1;
    power=param(Nparam); % Resouce decay rate
else
    power=0;
end
Nparam=Nparam+1;
kappa_r=param(Nparam); % Response variability
Nparam=Nparam+1;
K=param(Nparam); % Capacity
if ~isfield(Input,'Variants') % No Variants
    Input.Variants.Bias=0;
    Input.Variants.Swap=0;
    Input.Variants.BiasF=0;
    Input.Variants.PrecF=0;
end
if Input.Variants.Bias==0
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if Input.Variants.BiasF==0
    biasF=0; % Set bias as a consistent value
else
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
end
if Input.Variants.PrecF==0
    precF=0; % Set precision as consistent within each set size
else
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
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
if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
    sample_range=Data.sample_range;
    samples=Data.sample;
else
    samples=ones(1,length(errors));
    sample_range=1;
end
if Input.Variants.Swap==1
    errors_nt=Data.error_nt;
    if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
        samples_nt=Data.sample_nt;
    end
end
kappa_max=700; % Computational limit
SampleSeed=1000; % Monte Carlo seed
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input, SS_range); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
% LH function
if continuous==1
    
    p_error0=zeros(SampleSeed,length(errors));
    p_error0_NT=zeros(SampleSeed,length(errors));
    kappa=zeros(SampleSeed,length(SS_range),length(errors));
    kappa0_bar=exp(log(kappa1_bar)*(cosd(4*samples)).^precF);
    
    % MC Sampling
    for i_N=1:length(SS_range)
        kappa_bar=ones(SampleSeed,1)*kappa0_bar/(SS_range(i_N)).^power;
        kappa0=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
        kappa(:,i_N,:)=min(kappa0, kappa_max); % Constricted by the max kappa
    end
    
    for i_error=1:length(errors)
        N=SS(i_error);
        bias_cur=bias+biasF*cosd(4*samples(i_error)-90);
        error0=errors(i_error)+bias_cur; % Errors with bias
        % Von Mises distribution convoluted by kappa_r
        conv_kappa=sqrt(kappa(:,SS_range==N,i_error).^2+kappa_r^2+2*kappa(:,SS_range==N,i_error)*kappa_r.*cosd(error0));
        if K<N
            p_error0(:,i_error)=(1-K/N)*1/(error_range(2)-error_range(1))+...
                (K/N)*besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa(:,SS_range==N,i_error))*besseli0_fast(kappa_r));
        else
            p_error0(:,i_error)=besseli0_fast(conv_kappa)./(2*pi*besseli0_fast(kappa(:,SS_range==N,i_error))*besseli0_fast(kappa_r));
        end
        
        if Input.Variants.Swap==1
            if N==1
                p_error0_NT(:,i_error)=0;
            else
                p_temp_NT=zeros(SampleSeed,1);
                if K<N
                    for i_nt=1:N-1
                        error0_nt=errors_nt(i_nt, i_error)+bias_cur; % Errors with bias
                        conv_kappa_nt=sqrt(kappa(:,SS_range==N,i_error).^2+kappa_r^2+2*kappa(:,SS_range==N,i_error)*kappa_r.*cosd(error0_nt));
                        p_temp_NT=p_temp_NT+(1-K/N)*1/(error_range(2)-error_range(1))+...
                            (K/N)*besseli0_fast(conv_kappa_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N,i_error))*besseli0_fast(kappa_r));
                    end
                else
                    for i_nt=1:N-1
                        error0_nt=errors_nt(i_nt, i_error)+bias_cur; % Errors with bias
                        conv_kappa_nt=sqrt(kappa(:,SS_range==N,i_error).^2+kappa_r^2+2*kappa(:,SS_range==N,i_error)*kappa_r.*cosd(error0_nt));
                        p_temp_NT=p_temp_NT+besseli0_fast(conv_kappa_nt)./(2*pi*besseli0_fast(kappa(:,SS_range==N,i_error))*besseli0_fast(kappa_r));
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
            p_error0=p_error0(~isinf(sum(p_error0,2)),:,:);
            p_error(i_N,:,:)=mean(p_error0(~isnan(sum(p_error0,2)),:,:),1); % Find average across samples
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
            p_error0=p_error0(~isinf(sum(p_error0,2)),:,:);
            p_error(i_N,:,:)=mean(p_error0(~isnan(sum(p_error0,2)),:,:),1); % Find average across samples
            
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
        if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
            p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), sample_range==samples(i));
        else
            p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), 1);
        end
    end
    if Input.Variants.Swap==1
        for i=1:length(errors_nt)
            if SS(i)==1
                p_NT(i)=0;
            else
                for j=1:SS(i)-1
                    if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
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
    Output=log(LH);
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
function p=prior(param, Input, SS_range)

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
Nparam=Nparam+1;
kappa_r=param(Nparam); % Response variability
% Gamma prior for response noise
p0(Nparam)=gampdf(kappa_r,3,5);
Nparam=Nparam+1;
K=param(Nparam); % Capacity
% weibull prior for capacity
p0(Nparam)=wblpdf(K,3.5,3); % given that K is often 3~4
if ~isfield(Input,'Variants') % No Variants
    Input.Variants.Bias=0;
    Input.Variants.Swap=0;
    Input.Variants.BiasF=0;
    Input.Variants.PrecF=0;
end
if Input.Variants.Bias==1
    bias=param(Nparam); % Mean bias
    % Gaussian prior for bias
    p0(Nparam)=normpdf(bias, 0, 1);
end
if Input.Variants.BiasF==1
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
    % Gaussian prior for the fluctuation of bias
    p0(Nparam)=normpdf(biasF, 0, 5);
end
if Input.Variants.PrecF==1
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
    % Gaussian prior for the fluctuation of precision
    p0(Nparam)=normpdf(precF, 0, 1);
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
