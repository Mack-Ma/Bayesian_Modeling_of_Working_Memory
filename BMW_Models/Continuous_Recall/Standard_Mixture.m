%% Standard_Mixture
%
% Define the Standard Mixture model
% ------------
% Output=Standard_Mixture(param, Data, Input)
%
% ## Theory ##
% This model assumed that the fidelity of representation is unique in
% each set size & there's a fixed, discrete, item-based capacity limit.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% K, kappa_SS(:), (bias, biasF, precF, s)
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
% - Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations in visual working memory.
% Nature, 453(7192), 233.
% - van den Berg, R., Shin, H., Chou, W. C., George, R., & Ma, W. J. (2012).
% "Variability in encoding precision accounts for visual short-term memory limitations".
% Proceedings of the National Academy of Sciences, 109(22), 8780-8785.
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

function Output=Standard_Mixture(param, Data, Input)

if isfield(Input,'Variants') && any(strcmp(Input.Variants,'Category (Within-Item)'))
    Output=Categorical_Standard_Mixture_WI(param,Data,Input);
elseif isfield(Input,'Variants') && any(strcmp(Input.Variants,'Category (Between-Item)'))
    Output=Categorical_Standard_Mixture_BI(param,Data,Input);
else
    % Specify parameters
    K=param(1); % Capacity
    SS=Data.SS;
    SS_range=unique(SS);
    kappa_SS=param(2:1+length(SS_range)); % Precision
    Nparam=length(SS_range)+1;
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
    SS=Data.SS;
    SS_range=unique(SS);
    if ~any(strcmp(Input.Variants,'ContinuousK'))
        K=floor(K); % Note that capacity is a fixed, discrete value here
    end
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
    if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
        Prior=prior(param, Input, length(kappa_SS)); % get prior
    elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
        Prior=1; % uniform prior
    end
    
    if ~strcmp(Input.Output,'Prior')
        % LH function
        if continuous==1
            p_error=zeros(1,length(errors));
            p_error_NT=zeros(1,length(errors));
            kappa0=exp(log(kappa_SS')*(cosd(4*samples_cos)).^precF); % Flucutuative precision
            for i_error=1:length(errors)
                N=SS(i_error);
                error0=errors(i_error)+bias+biasF*cosd(4*samples_cos(i_error)-90);
                i_N=SS_range==N;
                if any(strcmp(Input.Variants,'ResponseNoise'))
                    if K<N % beyond capacity
                        conv_kappa=sqrt(kappa0(i_N,i_error).^2+kappa_r^2+2*kappa0(i_N,i_error)*kappa_r.*cosd(error0/(period/360)));
                        p_error(i_error)=besseli(0,conv_kappa)./(2*pi*besseli(0,kappa0(i_N,i_error))*besseli(0,kappa_r));
                    else % within capacity
                        p_error(i_error)=exp(kappa0(i_N,i_error).*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa0(i_N,i_error)));
                    end
                else
                    if K<N % beyond capacity
                        p_error(i_error)=(1-K/N)*1/(error_range(2)-error_range(1))+(K/N)*exp(kappa0(i_N,i_error).*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa0(i_N,i_error)));
                    else % within capacity
                        p_error(i_error)=exp(kappa0(i_N,i_error).*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa0(i_N,i_error)));
                    end
                end
                
                if any(strcmp(Input.Variants,'Swap'))
                    if N==1
                        p_error_NT(i_error)=0;
                    else
                        p_temp_NT=0;
                        if any(strcmp(Input.Variants,'ResponseNoise'))
                            for i_nt=1:N-1
                                error0_nt=errors_nt(i_error,i_nt)+bias+biasF*cosd(4*samples_nt_cos(i_error,i_nt)-90); % Errors with bias
                                conv_kappa_NT=sqrt(kappa0(i_N,i_error).^2+kappa_r^2+2*kappa0(i_N,i_error)*kappa_r.*cosd(error0_nt/(period/360)));
                                if K<N
                                    p_temp_NT=p_temp_NT+(1-K/N)*1/(error_range(2)-error_range(1))+(K/N)*besseli(0,conv_kappa_NT)./(2*pi*besseli(0,kappa0(i_N,i_error))*besseli(0,kappa_r));
                                else
                                    p_temp_NT=p_temp_NT+besseli(0,conv_kappa_NT)./(2*pi*besseli(0,kappa0(i_N,i_error))*besseli(0,kappa_r));
                                end
                            end
                        else
                            for i_nt=1:N-1
                                error0_nt=errors_nt(i_error,i_nt)+bias+biasF*cosd(4*samples_nt_cos(i_error,i_nt)-90); % Errors with bias
                                if K<N % beyond capacity
                                    p_temp_NT=p_temp_NT+(1-K/N)*1/(error_range(2)-error_range(1))+(K/N)*exp(kappa0(i_N,i_error).*cosd(error0_nt/(period/360)))./(2*pi*besseli(0,kappa0(i_N,i_error)));
                                else % within capacity
                                    p_temp_NT=p_temp_NT+exp(kappa0(i_N,i_error).*cosd(error0_nt/(period/360)))./(2*pi*besseli(0,kappa0(i_N,i_error)));
                                end
                            end
                        end
                        p_error_NT(i_error)=p_temp_NT/(N-1);
                    end
                end
            end
            p_T=(1-s)*p_error;
            p_NT=s*p_error_NT;
            p_LH=p_T+p_NT;
        else
            kappa0=zeros(length(SS_range),length(sample_range));
            for i_N=1:length(SS_range)
                kappa0(i_N,:)=exp(log(kappa_SS(i_N))*(cosd(4*sample_range)).^precF); % Flucutuative precision
            end
            bias_cur=bias+biasF*cosd(4*sample_range-90); % Fluctuative bias
            p_error=zeros(length(SS_range), length(error_range), length(sample_range));
            for i_N=1:length(SS_range)
                N=SS_range(i_N);
                if K<N % beyond capacity
                    if any(strcmp(Input.Variants,'ResponseNoise'))
                        for i_error=1:length(error_range)
                            error0=error_range(i_error)+bias_cur;
                            % Von Mises distribution convoluted by kappa_r
                            conv_kappa=sqrt(kappa0(i_N,:).^2+kappa_r^2+2*kappa0(i_N,:)*kappa_r.*cosd(error0/(period/360)));
                            p_error(i_N,i_error,:)=(1-K/N)*1/(2*period/360)/pi+...
                                (K/N)*besseli(0,conv_kappa)./(2*pi*besseli(0,kappa0(i_N,:))*besseli(0,kappa_r));
                        end
                    else
                        for i_error=1:length(error_range)
                            error0=error_range(i_error)+bias_cur;
                            p_error(i_N,i_error,:)=(1-K/N)*1/(2*period/360)/pi+...
                                (K/N)*exp(kappa0(i_N,:).*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa0(i_N,:)));
                        end
                    end
                    for i_sample=1:length(sample_range)
                        p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
                    end
                else
                    if any(strcmp(Input.Variants,'ResponseNoise'))
                        for i_error=1:length(error_range)
                            error0=error_range(i_error)+bias_cur;
                            % Von Mises distribution convoluted by kappa_r
                            conv_kappa=sqrt(kappa0(i_N,:).^2+kappa_r^2+2*kappa0(i_N,:)*kappa_r.*cosd(error0/(period/360)));
                            p_error(i_N,i_error,:)=besseli(0,conv_kappa)./(2*pi*besseli(0,kappa0(i_N,:))*besseli(0,kappa_r));
                        end
                    else
                        for i_error=1:length(error_range)
                            error0=error_range(i_error)+bias_cur;
                            p_error(i_N,i_error,:)=exp(kappa0(i_N,:).*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa0(i_N,:)));
                        end
                    end
                    for i_sample=1:length(sample_range)
                        p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample)); % Normalization
                    end
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
        
        if K>max(SS) % K cannot be larger than the set size
            LLH=realmax('double');
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
    
    if ~isstruct(Output) && (any(abs(Output))==Inf || any(isnan(LLH)))
        Output=realmax('double'); % LP should be a real value
    end
end

end

% Define prior
function p=prior(param, Input, Nkappa)

% Specify parameters
K=param(1); % Capacity
% weibull prior for capacity
p0(1)=wblpdf(K,3.5,3); % given that K is often 3~4
kappa_SS=param(2:1+Nkappa); % precision
% Gamma prior for precision
p0(2:1+Nkappa)=gampdf(kappa_SS,3,10);
Nparam=Nkappa+1;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'ResponseNoise'))
    Nparam=Nparam+1;
    kappa_r=param(Nparam);
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
