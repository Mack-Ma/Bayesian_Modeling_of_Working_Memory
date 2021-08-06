%% Target Confusability Competition
%
%
% Define the likelihood function for the TCC model
% ------------
% Output=TCC(param, Data, Input)
%
% ## Theory ##
% This model assumed that there are a limited number of resouce slots and
% the fidelity of each representation depends on the number of slots being
% allocated to the corresponding item.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% K, kappa_1, kappa_r, (bias, biasF, precF, s)
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
% - Zhang, W., & Luck, S. J. (2008). "Discrete fixed-resolution representations in visual working memory".
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
% 7/27/2021
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function Output=TCC(param, Data, Input)
% Specify parameters
SS=Data.SS;
SS_range=unique(SS);
dprime=param(1:length(SS_range)); % memory strength
%     Nparam=1;
%     if ~isfield(Input,'Variants') % No Variants
%         Input.Variants={};
%     end
%     if any(strcmp(Input.Variants,'ResponseNoise'))
%         Nparam=Nparam+1;
%         kappa_r=param(3); % Response precision
%     end

% Configuration
conf=Data.ConfusionVector;
samples=Data.sample;
responses=Data.response;
error_range=Data.error_range;
period=max(error_range)-min(error_range)+(error_range(2)-error_range(1));
errors=CircDist_BMW('Diff',responses,samples,period);

if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end


if ~strcmp(Input.Output,'Prior')
        
    p_error=zeros(length(SS_range), length(error_range));
    for s=1:length(SS_range)
        dp=dprime(s);
        mean_strength=dp.*conf;
        for i=1:length(error_range)
            p_error(s,i)=integral(@(x)strength_pdf(x,mean_strength(i),mean_strength(1:length(error_range)~=i)),-10,mean_strength(i))+...
                integral(@(x)strength_pdf(x,mean_strength(i),mean_strength(1:length(error_range)~=i)),mean_strength(i),10); % marginalize strength
        end
    end
    
    p_LH=zeros(1,length(errors));
    for i=1:length(errors)
        p_LH(i)=p_error(SS_range==SS(i),error_range==errors(i));
    end
    
    % LLH
    if isfield(Input,'PDF') && Input.PDF==1
        LLH.error=p_error; % PDF
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

function p=strength_pdf(s, mu1, mu2)
p=normpdf(s, mu1, 1); % pick one strength
for i=1:length(mu2)
    p_max=normcdf(s, mu2(i), 1); % probability of being larger than the other options
    p=p.*p_max;
end
end