%% Categorical Slots-plus-Averaging (Between-Variant)
%
% Log likelihood function of the Categorical Slots-plus-Averaging (Between-Variant) model
% ------------
% ## Theory ##
% This model assumed that there are a limited number of resouce slots and
% the fidelity of each representation depends on the number of slots being
% allocated to the corresponding item.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% - param
% K, kappa_1, kappa_r, kappa_c, p_c (bias, biasF, s)
% - Data
% Data.error (response-sample), Data.SS (set size), Data.error_c, (category-sample)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt, Data.error_nt_c)
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
% - Hardman, K. O., Vergauwe, E., & Ricker, T. J. (2017). Categorical working memory representations 
% are used in delayed estimation of continuous colors. Journal of Experimental Psychology: 
% Human Perception and Performance, 43(1), 30.
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


function LLH=Categorical_Slots_plus_Averaging_BV(param, Data, Input)

% Specify parameters
K=param(1);
kappa_1=param(2);
kappa_r=param(3);
kappa_c=param(4);
p_c=param(5);
Nparam=5;
if Input.Derivatives.Bias==0
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if Input.Derivatives.BiasF==1
    warning('biasF is not identified in the categorical models.')
end
if Input.Derivatives.PrecF==1
    warning('precF is not identified in the categorical models.')
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
SS=Data.SS;
SS_range=unique(SS);
errors_c=Data.error_c;
K=floor(K); % Note that capacity is a fixed, discrete value here
if Input.Derivatives.Swap==1
    errors_nt=Data.error_nt;
    errors_nt_c=Data.error_nt_c;
end

% LH function
p_error=zeros(length(SS_range), length(error_range));
p_error_c=zeros(length(SS_range), length(error_range));
for i_N=1:length(SS_range)
    N=SS_range(i_N);
    if K<N % Beyond capacity
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias;
            conv_1=sqrt((kappa_1).^2+(kappa_r)^2+2*kappa_1*kappa_r.*cosd(error0));
            conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0));
            p_error(i_N,i_error,1)=(1-K/N)*1/length(error_range)+...
                (K/N)*(besseli0_fast(conv_1)./(2*pi*besseli0_fast(kappa_1)*besseli0_fast(kappa_r)));
            p_error_c(i_N,i_error,1)=(1-K/N)*1/length(error_range)+...
                (K/N)*(besseli0_fast(conv_1_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r)));
        end
        p_error(i_N,:,1)=p_error(i_N,:,1)/sum(p_error(i_N,:,1));
        p_error_c(i_N,:,1)=p_error_c(i_N,:,1)/sum(p_error_c(i_N,:,1));
    else
        p_low=1-mod(K,N)/N;
        p_high=mod(K,N)/N;
        kappa_low=kappa_1*(floor(K/N)+1);
        kappa_high=kappa_1*floor(K/N);
        for i_error=1:length(error_range)
            error0=error_range(i_error)+bias;
            conv_1_c=sqrt((kappa_c).^2+(kappa_r)^2+2*kappa_c*kappa_r.*cosd(error0));
            conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0));
            conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0));
            p_error(i_N,i_error,:)=p_low*(besseli0_fast(conv_low)./(2*pi*besseli0_fast(kappa_low)*besseli0_fast(kappa_r)))+...
                p_high*(besseli0_fast(conv_high)./(2*pi*besseli0_fast(kappa_high)*besseli0_fast(kappa_r)));
            p_error_c(i_N,i_error,1)=(1-K/N)*1/length(error_range)+...
                (K/N)*(besseli0_fast(conv_1_c)./(2*pi*besseli0_fast(kappa_c)*besseli0_fast(kappa_r)));
        end
        p_error(i_N,:,1)=p_error(i_N,:,1)/sum(p_error(i_N,:,1));
        p_error_c(i_N,:,1)=p_error_c(i_N,:,1)/sum(p_error_c(i_N,:,1));
    end
end

% Calculate LH
p_T=zeros(1,length(errors));
p_NT=zeros(1,length(errors));
for i=1:length(errors)
    p_T(i)=(1-s)*((1-p_c)*p_error(SS_range==SS(i),error_range==errors(i),1)+...
        p_c*p_error_c(SS_range==SS(i),error_range==errors_c(i),1));
end
if Input.Derivatives.Swap==1
    for i=1:length(errors_nt)
        if SS(i)==1
            p_NT(i)=0;
        else
            for j=1:SS(i)-1
                p_NT(i)=p_NT(i)+s/(SS(i)-1)*((1-p_c)*p_error(SS_range==SS(i),error_range==errors_nt(j,i), 1)+...
                    p_c*p_error_c(SS_range==SS(i),error_range==errors_nt_c(j,i), 1));
            end
        end
    end
end
p_LH=p_T+p_NT; % Target + non-target LH

% LLH
LLH=-sum(log(p_LH));

if LLH==Inf || isnan(LLH)
    LLH=exp(666); % LLH should be a real value
end

end
