%% Category-Only
%
% Define the Category-Only model
% ------------
% Output=Category_Only(param, Data, Input)
%
% ## Theory ##
% This model assumed that all of the responses are made based on 
% the categorical working memory representations. The responses 
% follow Von Mises distributions that centered around the categories.
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% kappa_c (, kappa_r, bias, s)
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


function Output=Category_Only(param, Data, Input)

% Specify parameters
kappa_c=param(1); % Precision of categorical memory
Nparam=1;
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
if ~any(strcmp(Input.Variants,'Swap'))
    s=0; % No swap
    SS=ones(size(Data.sample));
else
    SS=Data.SS;
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
end

% Configuration
samples=Data.sample;
categories=Data.category;
error_range=Data.error_range;
if length(error_range)==2
    continuous=1;
    period=error_range(2)-error_range(1);
else
    continuous=0;
    period=max(error_range)-min(error_range)+error_range(2)-error_range(1);
end
errors_c=CircDist_BMW('Diff',categories,samples,period);
if any(strcmp(Input.Variants,'Swap'))
    samples_nt=Data.sample_nt;
    categories_nt=Data.category_nt;
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
        p_error=zeros(1,length(errors_c));
        p_error_NT=zeros(1,length(errors_c));
        for i_error=1:length(errors_c)
            error0_c=errors_c(i_error)+bias;
            if any(strcmp(Input.Variants,'ResponseNoise'))
                conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_c));
                p_error(i_error)=besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
            else
                p_error(i_error)=exp(kappa_c.*cosd(error0_c))./(2*pi*besseli(0,kappa_c));
            end
            
            if any(strcmp(Input.Variants,'Swap'))
                Nnt=SS(i_error)-1;
                if Nnt==0
                    p_error_NT(i_error)=0;
                else
                    p_temp_NT=0;
                    if any(strcmp(Input.Variants,'ResponseNoise'))
                        for i_nt=1:Nnt
                            error0_nt_c=errors_nt_c(i_error,i_nt)+bias;
                            conv_1_c_nt=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0_nt_c));
                            p_temp_NT=p_temp_NT+(besseli(0,conv_1_c_nt)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r)));
                        end
                    else
                        for i_nt=1:Nnt
                            error0_nt_c=errors_nt_c(i_error,i_nt)+bias;
                            p_temp_NT=p_temp_NT+(exp(kappa_c.*cosd(error0_nt_c))./(2*pi*besseli(0,kappa_c)));
                        end
                    end
                    p_error_NT(i_error)=p_temp_NT/Nnt;
                end
            end
        end
        p_T=(1-s)*p_error;
        p_NT=s*p_error_NT;
        p_LH=p_T+p_NT;
    else
        p_error_c=zeros(1, length(error_range));
        if any(strcmp(Input.Variants,'ResponseNoise'))
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias;
                conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0));
                p_error_c(i_error)=besseli(0,conv_1_c)./(2*pi*besseli(0,kappa_c)*besseli(0,kappa_r));
            end
        else
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias;
                p_error_c(i_error)=exp(kappa_c.*cosd(error0))./(2*pi*besseli(0,kappa_c));
            end
        end
        p_error_c=p_error_c/sum(p_error_c);
        
        % Calculate LH
        p_T=zeros(1,length(errors_c));
        p_NT=zeros(1,length(errors_c));
        for i=1:length(errors_c)
            p_T(i)=(1-s)*(p_error_c(error_range==errors_c(i)));
        end
        if any(strcmp(Input.Variants,'Swap'))
            for i=1:length(errors_nt_c)
                if SS(i)==1
                    p_NT(i)=0;
                else
                    for j=1:SS(i)-1
                        p_NT(i)=p_NT(i)+s/(SS(i)-1)*(p_error_c(error_range==errors_nt_c(i,j)));
                    end
                end
            end
        end
        p_LH=p_T+p_NT; % Target + non-target LH
    end
    
    % LLH
    if isfield(Input,'PDF') && Input.PDF==1
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
kappa_c=param(1); % categorical memory precision
% Gamma prior for categorical memory precision
p0(1)=gampdf(kappa_c,3,5);
Nparam=1;
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
% We think it's generally acceptable for prior definition.
p=1;
for i=1:Nparam
    p=p*p0(i);
end

end
