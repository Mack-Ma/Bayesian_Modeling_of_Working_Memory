%% Von Mises PDF
%
% Log likelihood function of Von Mises (circular normal) distribution
% ------------
% ## Input ##
% - param
% kappa, (bias)
% - Data
% Data.error (response-sample), Data.error_range
% - Input
% Input.Variants (options of Variants)
%
% ## Output ##
% Summed log Likelihood
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 12/8/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function LLH=VonMisesPDF(param, Data, Input)

% Specify parameters
kappa=param(1);
if isfield(Input.Variants, 'Bias') && Input.Variants.Bias==0
    bias=0; % Responses concentrate on samples
else
    bias=param(2);
end

% Configuration
errors=Data.error;
error_range=Data.error_range;

% Likelihood function
if ischar(error_range) && strcmp(error_range,'continuous') % When error is continuous
    p_error=zeros(1,length(errors));
    for i_error=1:length(errors)
        error0=bias+errors(i_error);
        p_error(i_error)=exp(kappa.*cosd(error0))./(2*pi*besseli0(kappa));
    end
    % Calculate LH
    p_LH=p_error;
elseif ~ischar(error_range) % When error is discrete 
    p_error=zeros(1,length(errors));
    for i_error=1:length(error_range)
        error0=bias+error_range(i_error);
        p_error(i_error)=exp(kappa.*cosd(error0))./(2*pi*besseli0(kappa));
    end
    % Normalization
    p_error=p_error/sum(p_error);
    % Calculate LH
    p_LH=zeros(1,length(errors));
    for i=1:length(errors)
        p_LH(i)=p_error(error_range==errors(i));
    end
else
    error('error_range is invalid...')
end

% LLH
LLH=-sum(log(p_LH)); % Negative LLH

if LLH==Inf || isnan(LLH)
    LLH=exp(666); % LLH should be a real value
end

end