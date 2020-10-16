%% Categorical Standard Mixture (Between-Item)
%
% Define the Categorical Standard Mixture (Between-Item) model
% ------------
% Output=Categorical_Standard_Mixture_BV(param, Data, Input)
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


function Output=Categorical_Standard_Mixture_BI(param, Data, Input)

% Specify parameters
K=param(1); % Capacity
SS=Data.SS;
SS_range=unique(SS);
kappa_SS=param(2:1+length(SS_range)); % Precision of continuous memory
Nparam=length(SS_range)+1;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
kappa_c=param(Nparam+1); % Precision of categorical memory
if any(strcmp(Input.Variants,'VariableCatWeight'))
    p_c=param((Nparam+2):(Nparam+1+length(SS_range)));
    Nparam=Nparam+1+length(SS_range);
else
    p_c=param(Nparam+2); % categorical weight
    Nparam=Nparam+2;
end
if any(strcmp(Input.Variants,'ResponseNoise'))
    Nparam=Nparam+1;
    kappa_r=param(Nparam);
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
SS=Data.SS;
SS_range=unique(SS);
if length(error_range)==2
    continuous=1;
    period=error_range(2)-error_range(1);
else
    continuous=0;
    period=max(error_range)-min(error_range)+error_range(2)-error_range(1);
end
errors=CircDist_BMW('Diff',responses,samples,period);
errors_c=CircDist_BMW('Diff',categories,samples,period);
if any(strcmp(Input.Variants,'Swap'))
    samples_nt=Data.sample_nt;
    categories_nt=Data.category_nt;
    errors_nt=CircDist_BMW('Diff',repmat(responses,[1,size(samples_nt,2)]),samples_nt,period);
    errors_nt_c=CircDist_BMW('Diff',categories_nt,samples_nt,period);
end
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
    % LH function
    if continuous==1
        p_error=zeros(1,length(errors));
        p_error_NT=zeros(1,length(errors));
        for i_error=1:length(errors)
            if any(strcmp(Input.Variants,'VariableCatWeight'))
                p_cc=p_c(SS_range==SS(i_error));
            else
                p_cc=p_c;
            end
            N=SS(i_error);
            error0=errors(i_error)+bias;
            error0_c=errors_c(i_error)+bias;
            i_N=SS_range==N;
            conv_SS=sqrt((kappa_SS).^2+(kappa_r)^2+2*kappa_SS*kappa_r.*cosd(error0));
            conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_c));
            if any(strcmp(Input.Variants,'ResponseNoise'))
                if K<N % Beyond capacity
                    p_error(i_error)=(1-K/N)*1/(error_range(2)-error_range(1))+...
                        (K/N)*((1-p_cc)*(besseli(0,conv_SS(i_N))./(2*pi*besseli(0,kappa_SS(i_N))*besseli(0,kappa_r)))+...
                        p_cc*(besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r))));
                else
                    p_error(i_error)=(1-p_cc)*(besseli(0,conv_SS(i_N))./(2*pi*besseli(0,kappa_SS(i_N))*besseli(0,kappa_r)))+...
                        p_cc*(besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                end
            else
                if K<N
                    p_error(i_error)=(1-K/N)*1/(error_range(2)-error_range(1))+...
                        (K/N)*((1-p_cc)*(exp(kappa_SS(i_N).*cosd(error0))./(2*pi*besseli(0,kappa_SS(i_N))))+...
                        p_cc*(exp(kappa_c.*cosd(error0_c))./(2*pi*besseli(0,kappa_c))));                        
                else
                    p_error(i_error)=(1-p_cc)*(exp(kappa_SS(i_N).*cosd(error0))./(2*pi*besseli(0,kappa_SS(i_N))))+...
                        p_cc*(exp(kappa_c.*cosd(error0_c))./(2*pi*besseli(0,kappa_c)));
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
                            conv_1_c_nt=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_nt_c));
                            conv_SS_nt=sqrt((kappa_SS).^2+(kappa_r)^2+2*kappa_SS*kappa_r.*cosd(error0_nt));
                            if K<N % beyond capacity
                                p_temp_NT=p_temp_NT+(1-K/N)*1/(error_range(2)-error_range(1))+...
                                    (K/N)*((1-p_cc)*(besseli(0,conv_SS_nt(i_N))./(2*pi*besseli(0,kappa_SS(i_N))*besseli(0,kappa_r))))+...
                                    p_cc*(besseli(0,conv_1_c_nt)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                            else % within capacity
                                p_temp_NT=p_temp_NT+(1-p_cc)*(besseli(0,conv_SS_nt(i_N))./(2*pi*besseli(0,kappa_SS(i_N))*besseli(0,kappa_r)))+...
                                    p_cc*(besseli(0,conv_1_c_nt)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                            end
                        end
                    else
                        for i_nt=1:N-1
                            error0_nt=errors_nt(i_error,i_nt)+bias; % Errors with bias
                            error0_nt_c=errors_nt_c(i_error,i_nt)+bias;
                            if K<N
                                p_temp_NT=p_temp_NT+(1-K/N)/(error_range(2)-error_range(1))+...
                                    (K/N)*((1-p_cc)*(exp(kappa_SS(i_N).*cosd(error0_nt))./(2*pi*besseli(0,kappa_SS(i_N))))+...
                                    p_cc*(exp(kappa_c.*cosd(error0_nt_c))./(2*pi*besseli(0,kappa_c))));
                            else
                                p_temp_NT=p_temp_NT+(1-p_cc)*(exp(kappa_SS(i_N).*cosd(error0_nt))./(2*pi*besseli(0,kappa_SS(i_N))))+...
                                    p_cc*(exp(kappa_c.*cosd(error0_nt_c))./(2*pi*besseli(0,kappa_c)));
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
    else % discrete mode
        p_error=zeros(length(SS_range), length(error_range));
        p_error_c=zeros(length(SS_range), length(error_range));
        for i_N=1:length(SS_range)
            N=SS_range(i_N);
            if any(strcmp(Input.Variants,'ResponseNoise'))
                if K<N % Beyond capacity
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        conv_SS=sqrt((kappa_SS).^2+(kappa_r)^2+2*kappa_SS*kappa_r.*cosd(error0));
                        conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0));
                        p_error(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                            (K/N)*(besseli(0,conv_SS(i_N))./(2*pi*besseli(0,kappa_SS(i_N))*besseli(0,kappa_r)));
                        p_error_c(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                            (K/N)*(besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                    end
                else
                    conv_SS=sqrt((kappa_SS).^2+(kappa_r)^2+2*kappa_SS*kappa_r.*cosd(error0));
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        p_error(i_N,i_error)=besseli(0,conv_SS(i_N))./(2*pi*besseli(0,kappa_SS(i_N))*besseli(0,kappa_r));
                        p_error_c(i_N,i_error)=besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
                    end
                end
            else
                if K<N
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        p_error(i_N,i_error)=(1-K/N)/length(error_range)+...
                            (K/N)*(exp(kappa_SS(i_N).*cosd(error0))./(2*pi*besseli(0,kappa_SS(i_N))));
                        p_error_c(i_N,i_error)=(1-K/N)/length(error_range)+...
                            (K/N)*(exp(kappa_c.*cosd(error0))./(2*pi*besseli(0,kappa_c)));
                    end
                else
                    for i_error=1:length(error_range)
                        error0=error_range(i_error)+bias;
                        p_error(i_N,i_error)=exp(kappa_SS(i_N).*cosd(error0))./(2*pi*besseli(0,kappa_SS(i_N)));
                        p_error_c(i_N,i_error)=exp(kappa_c.*cosd(error0))./(2*pi*besseli(0,kappa_c));
                    end
                end
            end
        end
        p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
        p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
        % Calculate LH
        p_T=zeros(1,length(errors));
        p_NT=zeros(1,length(errors));
        for i=1:length(errors)
            if any(strcmp(Input.Variants,'VariableCatWeight'))
                p_cc=p_c(SS_range==SS(i));
            else
                p_cc=p_c;
            end
            p_T(i)=(1-s)*((1-p_cc)*p_error(SS_range==SS(i),error_range==errors(i),1)+...
                p_cc*p_error_c(SS_range==SS(i),error_range==errors_c(i),1));
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
                        p_NT(i)=p_NT(i)+s/(SS(i)-1)*((1-p_c)*p_error(SS_range==SS(i),error_range==errors_nt(i,j), 1)+...
                            p_c*p_error_c(SS_range==SS(i),error_range==errors_nt_c(i,j), 1));
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
K=param(1); % Capacity
% weibull prior for capacity
p0(1)=wblpdf(K,3.5,3); % given that K is ofter 3~4
kappa_SS=param(2:1+Nkappa); % precision
% Gamma prior for precision
p0(2:1+Nkappa)=gampdf(kappa_SS,3,10);
Nparam=Nkappa+1;
kappa_c=param(Nparam+1); % categorical memory precision
% Gamma prior for categorical memory precision
p0(Nparam+1)=gampdf(kappa_c,3,5);
p_c=param(Nparam+2); % categorical weight
% Gaussian prior for the rate of categorical memory
p0(Nparam+2)=normpdf(p_c, 0.5, 1);
Nparam=Nparam+2;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
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
