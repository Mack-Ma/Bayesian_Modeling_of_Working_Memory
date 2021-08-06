%% Categorical Slots-plus-Averaging (Between-Item)
%
% Define the Categorical Slots-plus-Averaging (Between-Item) model
% ------------
% Output=Categorical_Slots_plus_Averaging_BI(param, Data, Input)
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


function Output=Categorical_Slots_plus_Averaging_BI(param, Data, Input)

% Specify parameters
K=param(1); % Capacity
kappa_1=param(2); % Unit precision
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
SS=Data.SS;
SS_range=unique(SS);
kappa_c=param(3); % Precision of categorical memory
if any(strcmp(Input.Variants,'VariableCatWeight'))
    p_c=param(4:4+length(SS_range)-1);
    Nparam=3+length(SS_range);
else
    p_c=param(4); % Weight of categorical memory
    Nparam=4;
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
% if Input.Variants.BiasF==1
%     warning('biasF is not identified in the categorical models.')
% end
% if Input.Variants.PrecF==1
%     warning('precF is not identified in the categorical models.')
% end
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
    continuous=1;
    period=error_range(2)-error_range(1);
else
    continuous=0;
    period=max(error_range)-min(error_range)+error_range(2)-error_range(1);
end
errors=CircDist_BMW('Diff',responses,samples,period);
errors_c=CircDist_BMW('Diff',responses,categories,period);
if any(strcmp(Input.Variants,'Swap'))
    samples_nt=Data.sample_nt;
    categories_nt=Data.category_nt;
    errors_nt=CircDist_BMW('Diff',repmat(responses,[1,size(samples_nt,2)]),samples_nt,period);
    errors_nt_c=CircDist_BMW('Diff',repmat(responses,[1,size(samples_nt,2)]),categories_nt,period);
end
if ~any(strcmp(Input.Variants,'ContinuousK'))
    K=floor(K); % Note that capacity is a fixed, discrete value here
end
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Data, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
    % LH function
    if continuous==1
        p_error=zeros(1,length(errors));
        p_error_NT=zeros(1,length(errors));
        for i_error=1:length(errors)
            N=SS(i_error);
            if any(strcmp(Input.Variants,'VariableCatWeight'))
                p_cc=p_c(SS_range==N);
            else
                p_cc=p_c;
            end
            error0=errors(i_error)+bias;
            error0_c=errors_c(i_error)+bias;
            p_low=1-mod(K,N)/N;
            p_high=mod(K,N)/N;
            kappa_high=kappa_1*(floor(K/N)+1);
            kappa_low=kappa_1*floor(K/N);
            if any(strcmp(Input.Variants,'ResponseNoise'))
                conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_c/(period/360)));
                if K<N % Beyond capacity
                    conv_1=sqrt((kappa_1).^2+(kappa_r)^2+2*kappa_1*kappa_r.*cosd(error0/(period/360)));
                    p_error(i_error)=(1-K/N)*1/(2*period/360)/pi+...
                        (K/N)*((1-p_cc)*(besseli(0,conv_1)./(2*pi*besseli(0,kappa_1)*besseli(0,kappa_r)))+...
                        p_cc*(besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r))));
                else
                    conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0/(period/360)));
                    conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0/(period/360)));
                    p_error(i_error)=(1-p_cc)*(p_low*(besseli(0,conv_low)./(2*pi*besseli(0,kappa_low)*besseli(0,kappa_r)))+...
                        p_high*(besseli(0,conv_high)./(2*pi*besseli(0,kappa_high)*besseli(0,kappa_r))))+...
                        p_cc*besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
                end
            else
                if K<N
                    p_error(i_error)=(1-K/N)*1/(2*period/360)/pi+...
                        (K/N)*((1-p_cc)*(exp(kappa_1.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_1)))+...
                        p_cc*(exp(kappa_c.*cosd(error0_c/(period/360)))./(2*pi*besseli(0,kappa_c))));
                else
                    p_error(i_error)=(1-p_cc)*(p_low*(exp(kappa_low.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_low)))+...
                        p_high*(exp(kappa_high.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_high))))+...
                        p_cc*(exp(kappa_c.*cosd(error0_c/(period/360)))./(2*pi*besseli(0,kappa_c)));
                end
            end
            
            if any(strcmp(Input.Variants,'Swap'))
                if N==1
                    p_error_NT(i_error)=0;
                else
                    p_temp_NT=0;
                    if any(strcmp(Input.Variants,'ResponseNoise'))
                        for i_nt=1:N-1
                            error0_nt=errors_nt(i_error,i_nt)+bias; % Errors with bias
                            error0_nt_c=errors_nt_c(i_error,i_nt)+bias;
                            conv_1_c_nt=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_nt_c/(period/360)));
                            if K<N % beyond capacity
                                conv_1=sqrt((kappa_1).^2+(kappa_r)^2+2*kappa_1*kappa_r.*cosd(error0_nt/(period/360)));
                                p_temp_NT=p_temp_NT+(1-K/N)*1/(2*period/360)/pi+...
                                    (K/N)*((1-p_cc)*(besseli(0,conv_1)./(2*pi*besseli(0,kappa_1)*besseli(0,kappa_r))))+...
                                    p_cc*(besseli(0,conv_1_c_nt)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                            else % within capacity
                                conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0_nt/(period/360)));
                                conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0_nt/(period/360)));
                                p_temp_NT=p_temp_NT+(1-p_cc)*(p_low*(besseli(0,conv_low)./(2*pi*besseli(0,kappa_low)*besseli(0,kappa_r)))+...
                                    p_high*(besseli(0,conv_high)./(2*pi*besseli(0,kappa_high)*besseli(0,kappa_r))))+...
                                    p_cc*(besseli(0,conv_1_c_nt)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                            end
                        end
                    else
                        for i_nt=1:N-1
                            error0_nt=errors_nt(i_error,i_nt)+bias; % Errors with bias
                            error0_nt_c=errors_nt_c(i_error,i_nt)+bias;
                            if K<N
                                p_temp_NT=p_temp_NT+(1-K/N)*1/(2*period/360)/pi+...
                                    (K/N)*((1-p_cc)*(exp(kappa_.*cosd(error0_nt/(period/360)))./(2*pi*besseli(0,kappa_c)))+...
                                    p_cc*(exp(kappa_c.*cosd(error0_nt_c/(period/360)))./(2*pi*besseli(0,kappa_c))));
                            else
                                p_temp_NT=p_temp_NT+(1-p_cc)*(p_low*(exp(kappa_low.*cosd(error0_nt/(period/360)))./(2*pi*besseli(0,kappa_low)))+...
                                    p_high*(exp(kappa_high.*cosd(error0_nt/(period/360)))./(2*pi*besseli(0,kappa_high))))+...
                                    p_cc*(exp(kappa_c.*cosd(error0_nt_c/(period/360)))./(2*pi*besseli(0,kappa_c)));
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
        p_error=zeros(length(SS_range), length(error_range));
        p_error_c=zeros(length(SS_range), length(error_range));
        for i_N=1:length(SS_range)
            N=SS_range(i_N);
            if K<N % Beyond capacity
                if any(strcmp(Input.Variants,'ResponseNoise'))
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        conv_1=sqrt((kappa_1).^2+(kappa_r)^2+2*kappa_1*kappa_r.*cosd(error0/(period/360)));
                        conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0/(period/360)));
                        p_error(i_N,i_error)=besseli(0,conv_1)./(2*pi*besseli(0,kappa_1)*besseli(0,kappa_r));
                        p_error_c(i_N,i_error)=besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
                    end
                else
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        %                         p_error(i_N,i_error)=(1-K/N)/(2*period/360)/pi+...
                        %                             (K/N)*(exp(kappa_1.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_1)));
                        p_error(i_N,i_error)=exp(kappa_1.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_1));
                        %                         p_error_c(i_N,i_error)=(1-K/N)/(2*period/360)/pi+...
                        %                             (K/N)*(exp(kappa_c.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_c)));
                        p_error_c(i_N,i_error)=exp(kappa_c.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_c));
                    end
                end
                p_error(i_N,:)=(K/N)*p_error(i_N,:)/sum(p_error(i_N,:))+(1-K/N)/length(error_range);
                p_error_c(i_N,:)=(K/N)*p_error_c(i_N,:)/sum(p_error_c(i_N,:))+(1-K/N)/length(error_range);
            else
                p_low=1-mod(K,N)/N;
                p_high=mod(K,N)/N;
                kappa_high=kappa_1*(floor(K/N)+1);
                kappa_low=kappa_1*floor(K/N);
                if any(strcmp(Input.Variants,'ResponseNoise'))
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0/(period/360)));
                        conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0/(period/360)));
                        conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0/(period/360)));
                        p_error(i_N,i_error)=p_low*(besseli(0,conv_low)./(2*pi*besseli(0,kappa_low)*besseli(0,kappa_r)))+...
                            p_high*(besseli(0,conv_high)./(2*pi*besseli(0,kappa_high)*besseli(0,kappa_r)));
                        p_error_c(i_N,i_error)=besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
                    end
                else
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        p_error(i_N, i_error)=p_low*(exp(kappa_low.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_low)))+...
                            p_high*(exp(kappa_high.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_high)));
                        p_error_c(i_N, i_error)=exp(kappa_c.*cosd(error0/(period/360)))./(2*pi*besseli(0,kappa_c));
                    end
                end
                p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
                p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
            end
        end
        
        % Calculate LH
        p_T=zeros(1,length(errors));
        p_NT=zeros(1,length(errors));
        for i=1:length(errors)
            if any(strcmp(Input.Variants,'VariableCatWeight'))
                p_cc=p_c(SS_range==SS(i));
            else
                p_cc=p_c;
            end
            p_T(i)=(1-s)*((1-p_cc)*p_error(SS_range==SS(i),error_range==errors(i))+...
                p_cc*p_error_c(SS_range==SS(i),error_range==errors_c(i)));
        end
        if any(strcmp(Input.Variants,'Swap'))
            for i=1:length(errors_nt)
                if any(strcmp(Input.Variants,'VariableCatWeight'))
                    p_cc=p_c(SS_range==SS(i));
                else
                    p_cc=p_c;
                end
                if SS(i)==1
                    p_NT(i)=0;
                else
                    for j=1:SS(i)-1
                        p_NT(i)=p_NT(i)+s/(SS(i)-1)*((1-p_cc)*p_error(SS_range==SS(i),error_range==errors_nt(i,j), 1)+...
                            p_cc*p_error_c(SS_range==SS(i),error_range==errors_nt_c(i,j), 1));
                    end
                end
            end
        end
        p_LH=p_T+p_NT; % Target + non-target LH
    end
    
    % LLH
    if isfield(Input,'PDF') && Input.PDF==1
        LLH.error=p_error; % PDF of continuous error
        LLH.error_c=p_error_c; % PDF of categorical error
    else
        LLH=-sum(log(p_LH)); % Negative LLH
    end    
    
    % Posterior
    LP=-log(Prior)+LLH; % likelihood*prior
    
    if LP==Inf || isnan(LP)
        LP=realmax('double'); % Output should be a real value
    end
    
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

end

% Define prior
function p=prior(param, Data, Input)

% Specify parameters
K=param(1); % Capacity
% weibull prior for capacity
p0(1)=wblpdf(K,3.5,3); % given that K is ofter 3~4
kappa_1=param(2); % Unit resource
% Gamma prior for unit resource
p0(2)=gampdf(kappa_1,3,5);
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
SS=Data.SS;
SS_range=unique(SS);
if any(strcmp(Input.Variants,'VariableCatWeight'))
    kappa_c=param(3); % categorical memory precision
    % Gamma prior for categorical memory precision
    p0(3)=gampdf(kappa_c,3,5);
    p_c=param(4:4+length(SS_range)-1); % categorical weight
    % Gaussian prior for the rate of categorical memory
    p0(4:4+length(SS_range)-1)=normpdf(p_c, 0.5, 1);
    Nparam=3+length(SS_range);
else
    kappa_c=param(3); % categorical memory precision
    % Gamma prior for categorical memory precision
    p0(3)=gampdf(kappa_c,3,5);
    p_c=param(4); % categorical weight
    % Gaussian prior for the rate of categorical memory
    p0(4)=normpdf(p_c, 0.5, 1);
    Nparam=4;
end
if any(strcmp(Input.Variants,'ResponseNoise'))
    Nparam=Nparam+1;
    kappa_r=param(Nparam); % Response precision
    % Gamma prior for response noise
    p0(Nparam)=gampdf(kappa_r,3,5);
end
if any(strcmp(Input.Variants,'Bias'))
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
    % Gaussian prior for bias
    p0(Nparam)=normpdf(bias, 0, 1);
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
