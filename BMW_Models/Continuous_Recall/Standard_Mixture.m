%% Standard_Mixture
%
% Log likelihood function of the Standard Mixture model
% ------------
% ## Theory ##
% This model assumed that the fidelity of representations is difference in
% each set size & there's a fixed, discrete, item-based capacity limit.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% - param
% K, kappa_S(:), (bias, biasF, s)
% - Data
% Data.error (response-sample), Data.SS (set size)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt)
% - Input
% Input.Derivatives (options of derivatives)
%
% ## Output ##
% Summed log Likelihood
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

function LLH=Standard_Mixture(param, Data, Input)

% Specify parameters
K=param(1);
SS=Data.SS;
SS_range=unique(SS);
kappa_SS=param(2:1+length(SS_range));
Nparam=length(SS_range)+1;
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
K=floor(K); % Note that capacity is a fixed, discrete value here
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

% LH function
kappa0=zeros(length(SS_range),length(sample_range));
for i_N=1:length(SS_range)
    kappa0(i_N,:)=exp(log(kappa_SS(i_N))*(cosd(4*sample_range)).^precF); % Flucutuative precision
end
bias_cur=bias+biasF*cosd(4*sample_range-90); % Fluctuative bias
p_error=zeros(length(SS_range), length(error_range), length(sample_range));
for i_N=1:length(SS_range)
    N=SS_range(i_N);
    if K<N % beyond capacity
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur;
            p_error(i_N,i_error,:)=(1-K/N)*1/length(error_range)+...
                (K/N)*exp(kappa0(i_N,:).*cosd(error0))./(2*pi*besseli0_fast(kappa0(i_N,:)));
        end
        for i_sample=1:length(sample_range)
            p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
        end
    else
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias_cur;
            p_error(i_N,i_error,:)=exp(kappa0(i_N,:).*cosd(error0))./(2*pi*besseli0_fast(kappa0(i_N,:)));
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

% LLH
LLH=-sum(log(p_LH)); % Negative LLH

if LLH==Inf || isnan(LLH)
    LLH=exp(666); % LLH should be a real value
end

end

