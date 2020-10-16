%% Categorical Slots-plus-Averaging (Within-Item)
%
% Define the Categorical Slots-plus-Averaging (Within-Item) model
% ------------
% Output=Categorical_Slots_plus_Averaging_WI(param, Data, Input)
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


function Output=Categorical_Slots_plus_Averaging_WI(param, Data, Input)

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
    gamma_c=param(4:3+length(SS_range)); % Scaling the weight of categorical memory
    Nparam=3+length(SS_range);
else
    gamma_c=param(4); % Scaling the weight of categorical memory
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
    period=error_range(2)-error_range(1);
else
    period=max(error_range)-min(error_range)+(error_range(2)-error_range(1));
end
errors=CircDist_BMW('Diff',responses,samples,period);
errors_c=CircDist_BMW('Diff',responses,categories,period);
response_range=Data.response_range;
sample_range=Data.sample_range;
category_range=unique(Data.category);
if ~any(strcmp(Input.Variants,'ContinuousK'))
    K=floor(K); % Note that capacity is a fixed, discrete value here
end
if any(strcmp(Input.Variants,'Swap'))
    samples_nt=Data.sample_nt;
    categories_nt=Data.category_nt;
    errors_nt=CircDist_BMW('Diff',repmat(responses,1,size(samples_nt,2)),samples_nt,period);
    errors_nt_c=CircDist_BMW('Diff',repmat(responses,1,size(samples_nt,2)),categories_nt,period);
end
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Data, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

if ~strcmp(Input.Output,'Prior')
    % LH function
    p_error=zeros(length(SS_range), length(error_range));
    p_error_c=zeros(length(SS_range), length(error_range));
    for i_N=1:length(SS_range)
        N=SS_range(i_N);
        if K<N % Beyond capacity
            if any(strcmp(Input.Variants,'ResponseNoise'))
                for i_error=1:length(error_range)
                    error0=error_range(i_error)+bias;
                    conv_1=sqrt((kappa_1).^2+(kappa_r)^2+2*kappa_1*kappa_r.*cosd(error0));
                    conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0));
                    p_error(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                        (K/N)*(besseli(0,conv_1)./(2*pi*besseli(0,kappa_1)*besseli(0,kappa_r)));
                    p_error_c(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                        (K/N)*(besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                end
            else
                for i_error=1:length(error_range)
                    error0=error_range(i_error)+bias;
                    p_error(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                        (K/N)*(exp(kappa_1.*cosd(error0))./(2*pi*besseli(0,kappa_1)));
                    p_error_c(i_N,i_error)=(1-K/N)*1/length(error_range)+...
                        (K/N)*(exp(kappa_c.*cosd(error0))./(2*pi*besseli(0,kappa_c)));
                end
            end
            p_error(i_N,:)=p_error(i_N,:,1)/sum(p_error(i_N,:,1));
            p_error_c(i_N,:)=p_error_c(i_N,:,1)/sum(p_error_c(i_N,:,1));
        else
            p_low=1-mod(K,N)/N;
            p_high=mod(K,N)/N;
            kappa_low=kappa_1*(floor(K/N)+1);
            kappa_high=kappa_1*floor(K/N);
            if any(strcmp(Input.Variants,'ResponseNoise'))
                for i_error=1:length(error_range)
                    error0=error_range(i_error)+bias;
                    conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0));
                    conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0));
                    conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0));
                    p_error(i_N,i_error)=p_low*(besseli(0,conv_low)./(2*pi*besseli(0,kappa_low)*besseli(0,kappa_r)))+...
                        p_high*(besseli(0,conv_high)./(2*pi*besseli(0,kappa_high)*besseli(0,kappa_r)));
                    p_error_c(i_N,i_error)=besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
                end
            else
                for i_error=1:length(error_range)
                    error0=error_range(i_error)+bias;
                    p_error(i_N,i_error)=p_low*(exp(kappa_low.*cosd(error0))./(2*pi*besseli(0,kappa_low)))+...
                        p_high*(exp(kappa_high.*cosd(error0))./(2*pi*besseli(0,kappa_high)));
                    p_error_c(i_N,i_error)=exp(kappa_c.*cosd(error0))./(2*pi*besseli(0,kappa_c));
                end
            end
            p_error(i_N,:)=p_error(i_N,:)/sum(p_error(i_N,:));
            p_error_c(i_N,:)=p_error_c(i_N,:)/sum(p_error_c(i_N,:));
        end
    end
    
    % calculate the normalization constant
    nc=zeros(length(SS_range),length(sample_range));
    period_sample=sample_range(end)+sample_range(2)-2*sample_range(1);
    for i_N=1:length(SS_range)
        if any(strcmp(Input.Variants,'VariableCatWeight'))
            gamma_cc=gamma_c(i_N);
        else
            gamma_cc=gamma_c;
        end
        for i_s=1:length(sample_range)
            S=sample_range(i_s);
            % find the current category
            [~,C_ind]=min(abs(CircDist_BMW('Diff',repmat(S,size(category_range,1),1),category_range,period_sample)));
            C=category_range(C_ind);
            e_s=CircDist_BMW('Diff',response_range,ones(size(response_range))*S,period_sample); % continuous error
            e_c=CircDist_BMW('Diff',response_range,ones(size(response_range))*C,period_sample); % categorical error
            s_ind=interp1(error_range,1:length(error_range),e_s);
            c_ind=interp1(error_range,1:length(error_range),e_c);
            nc(i_N,i_s)=sum(p_error(i_N,s_ind).*(p_error_c(i_N,c_ind)).^gamma_cc);
        end
    end
    
    % Calculate LH
    p_T=zeros(1,length(errors));
    p_NT=zeros(1,length(errors));
    for i=1:length(errors)
        if any(strcmp(Input.Variants,'VariableCatWeight'))
            gamma_cc=gamma_c(SS_range==SS(i));
        else
            gamma_cc=gamma_c;
        end
        p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i))*(p_error_c(SS_range==SS(i),error_range==errors_c(i))).^gamma_cc/...
            nc(SS_range==SS(i),sample_range==samples(i));
    end
    if any(strcmp(Input.Variants,'Swap'))
        for i=1:length(errors_nt)
            if any(strcmp(Input.Variants,'VariableCatWeight'))
                gamma_cc=gamma_c(SS_range==SS(i));
            else
                gamma_cc=gamma_c;
            end
            if SS(i)==1
                p_NT(i)=0;
            else
                for j=1:SS(i)-1
                    p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(i,j)).*(p_error_c(SS_range==SS(i),error_range==errors_nt_c(i,j))).^gamma_cc./nc(SS_range==SS(i),sample_range==samples_nt(i,j));
                end
            end
        end
    end
    p_LH=p_T+p_NT; % Target + non-target LH
    
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
    gamma_c=param(4:3+length(SS_range)); % categorical weight
    % Gamma prior for the weight of categorical memory
    p0(4:3+length(SS_range))=gampdf(gamma_c, 2, 2);
    Nparam=3+length(SS_range);
else
    kappa_c=param(3); % categorical memory precision
    % Gamma prior for categorical memory precision
    p0(3)=gampdf(kappa_c,3,5);
    gamma_c=param(4); % categorical weight
    % Gamma prior for the weight of categorical memory
    p0(4)=gampdf(gamma_c, 2, 2);
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
