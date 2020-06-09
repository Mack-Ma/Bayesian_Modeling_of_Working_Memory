%% Variable Precision
%
% Define the Variable Precision model
% ------------
% Output=Variable_Precision(param, Data, Input)
%
% ## Theory ##
% This model assumed that there's no item-based capacity limit and
% the relationship between memory fidelity and set size follows a power-law
% function. In addition, memory resource is supposed to have 
% trial-by-trial, item-by-item variation.
% Response follows a Von Mises distribution in which the shape parameter
% follows the Gamma distribution
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% kappa1_bar, tau, (power,) kappa_r, (bias, biasF, precF, s, (kappa_c, p_c)/(kappa_c, eps_c))
%
% - Data
% Data.sample, Data.response, Data.SS (set size), Data.error_range
% (Data.sample_range, Data.sample_nt, Data.response_nt)
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
%       'Prior', only output the prior density
%       'LLH', output log likelihood
%       'LP', output log posterior density
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

function Output=Variable_Precision(param, Data, Input)

if isfield(Input,'Variants') && any(strcmp(Input.Variants,'Category (Within-Item)'))
    Output=Categorical_Variable_Precision_WI(param,Data,Input);
elseif isfield(Input,'Variants') && any(strcmp(Input.Variants,'Category (Between-Item)'))
    Output=Categorical_Variable_Precision_BI(param,Data,Input);
else
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
if ~any(strcmp(Input.Variants,'BiasF'))
    biasF=0; % Set bias as a consistent value
else
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
end
if ~any(strcmp(Input.Variants,'PrecF'))
    precF=0; % Set precision as consistent within each set size
else
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
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
error_range=Data.error_range;
if length(error_range)==2
    continuous=1;
    period=error_range(2)-error_range(1);
else
    continuous=0;
    period=max(error_range)-min(error_range)+(error_range(2)-error_range(1));
end
errors=CircDist_BMW('Diff',responses,samples,period);
if any(strcmp(Input.Variants,'BiasF')) || any(strcmp(Input.Variants,'PrecF'))
    if isfield(Data,'sample_range')
        sample_range=Data.sample_range;
    end
    samples_cos=samples;
else
    samples_cos=ones(1,length(errors));
    sample_range=1;
end
if any(strcmp(Input.Variants,'Swap'))
    samples_nt=Data.sample_nt;
    errors_nt=CircDist_BMW('Diff',repmat(responses,[1,size(samples_nt,2)]),samples_nt,period);
    if any(strcmp(Input.Variants,'BiasF')) || any(strcmp(Input.Variants,'PrecF'))
        samples_nt_cos=samples_nt;
    else
        samples_nt_cos=ones(size(samples_nt,1),size(samples_nt,2));
    end
end
kappa_max=700; % Computational limit
SampleSeed=2000; % Monte Carlo seed
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input, SS_range); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
    % LH function
    if continuous==1 % Continuous mode
        
        p_error0=zeros(SampleSeed,length(errors));
        p_error0_NT=zeros(SampleSeed,length(errors));
        kappa=zeros(SampleSeed,length(SS_range),length(errors));
        kappa0_bar=exp(log(kappa1_bar)*(cosd(4*samples_cos)).^precF);
        
        % MC Sampling
        for i_N=1:length(SS_range)
            kappa_bar=ones(SampleSeed,1)*kappa0_bar/(SS_range(i_N)).^power;
            kappa0=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
            kappa(:,i_N,:)=min(kappa0, kappa_max); % Constricted by the max kappa
        end
        
        for i_error=1:length(errors)
            N=SS(i_error);
            bias_cur=bias+biasF*cosd(4*samples_cos(i_error)-90);
            error0=errors(i_error)+bias_cur; % Errors with bias
            if any(strcmp(Input.Variants,'ResponseNoise'))
                % Von Mises distribution convoluted by kappa_r
                conv_kappa=sqrt(kappa(:,SS_range==N,i_error).^2+kappa_r^2+2*kappa(:,SS_range==N,i_error)*kappa_r.*cosd(error0));
                p_error0(:,i_error)=besseli(0,conv_kappa)./(2*pi*besseli(0,kappa(:,SS_range==N,i_error))*besseli(0,kappa_r));
            else
                p_error0(:,i_error)=exp(kappa(:,SS_range==N,i_error).*cosd(error0))./(2*pi*besseli(0,kappa(:,SS_range==N,i_error)));
            end
            if any(strcmp(Input.Variants,'Swap'))
                if N==1
                    p_error0_NT(:,i_error)=0;
                else
                    p_temp_NT=zeros(SampleSeed,1);
                    if any(strcmp(Input.Variants,'ResponseNoise'))
                        for i_nt=1:N-1
                            error0_nt=errors_nt(i_error,i_nt)+samples_nt_cos(i_error,i_nt); % Errors with bias
                            conv_kappa_nt=sqrt(kappa(:,SS_range==N,i_error).^2+kappa_r^2+2*kappa(:,SS_range==N,i_error)*kappa_r.*cosd(error0_nt));
                            p_temp_NT=p_temp_NT+besseli(0,conv_kappa_nt)./(2*pi*besseli(0,kappa(:,SS_range==N,i_error))*besseli(0,kappa_r));
                        end
                    else
                        for i_nt=1:N-1
                            error0_nt=errors_nt(i_error,i_nt)+samples_nt_cos(i_error,i_nt); % Errors with bias
                            p_temp_NT=p_temp_NT+exp(kappa(:,SS_range==N,i_error).*cosd(error0_nt))./(2*pi*besseli(0,kappa(:,SS_range==N,i_error)));
                        end
                    end
                    p_error0_NT(:,i_error)=p_temp_NT/(N-1);
                end
            end
        end
        p_error0=p_error0(~isinf(sum(p_error0,2)),:);
        p_T=(1-s)*mean(p_error0(~isnan(sum(p_error0,2)),:),1);
        p_error0_NT=p_error0_NT(~isinf(sum(p_error0_NT,2)),:);
        p_NT=s*mean(p_error0_NT(~isnan(sum(p_error0_NT,2)),:),1);
        p_LH=p_T+p_NT;
        
    else % Discrete mode
        
        kappa0_bar=exp(log(kappa1_bar)*(cosd(4*sample_range)).^precF);
        bias_cur=bias+biasF*cosd(4*sample_range-90);
        bias_cur=ones(SampleSeed,1)*bias_cur;
        p_error=zeros(length(SS_range),length(error_range), length(sample_range));
        for i_N=1:length(SS_range)
            N=SS_range(i_N);
            
            % MC Sampling
            kappa_bar=ones(SampleSeed,1)*kappa0_bar/(N).^power;
            kappa=gamrnd(kappa_bar/tau, tau); % Sample from gamma distribution
            kappa=min(kappa, kappa_max); % Constricted by the max kappa
            
            p_error0=zeros(SampleSeed,length(error_range),length(sample_range));
            if any(strcmp(Input.Variants,'ResponseNoise'))
                for i_error=1:length(error_range)
                    error0=error_range(i_error)+bias_cur; % Errors with bias
                    % Von Mises distribution convoluted by kappa_r
                    conv_kappa=sqrt(kappa.^2+kappa_r^2+2*kappa*kappa_r.*cosd(error0));
                    p_error0(:,i_error,:)=besseli(0,conv_kappa)./(2*pi*besseli(0,kappa)*besseli(0,kappa_r));
                end
            else
                for i_error=1:length(error_range)
                    error0=error_range(i_error)+bias_cur;
                    p_error0(:,i_error,:)=exp(kappa.*cosd(error0))./(2*pi*besseli(0,kappa));
                end
            end
            p_error0=p_error0(~isinf(sum(p_error0,2)),:,:);
            p_error(i_N,:,:)=mean(p_error0(~isnan(sum(p_error0,2)),:,:),1); % Find average across samples
            
            % Normalization
            for i_sample=1:length(sample_range)
                p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
            end
        end
        % Calculate LH
        p_T=zeros(1,length(errors));
        p_NT=zeros(1,length(errors));
        for i=1:length(errors)
            if any(strcmp(Input.Variants,'BiasF')) || any(strcmp(Input.Variants,'PrecF'))
                p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), sample_range==samples_cos(i));
            else
                p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), 1);
            end
        end
        if any(strcmp(Input.Variants,'Swap'))
            for i=1:length(errors_nt)
                if SS(i)==1
                    p_NT(i)=0;
                else
                    for j=1:SS(i)-1
                        if any(strcmp(Input.Variants,'BiasF')) || any(strcmp(Input.Variants,'PrecF'))
                            p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(i,j), sample_range==samples_nt_cos(i,j));
                        else
                            p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(i,j), 1);
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
    Output=log(p_LH);
elseif strcmp(Input.Output,'All')
    Output.LP=LP;
    Output.LLH=LLH;
    Output.Prior=Prior;
    Output.LPPD=log(p_LH);
end

if ~isstruct(Output) && (any(isinf(abs(Output))) || any(isnan(Output)))
    Output=realmax('double'); % Output should be a real value
end

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
