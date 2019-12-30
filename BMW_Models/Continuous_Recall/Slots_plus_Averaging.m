%% Slots-plus-Averaging
%
% Log likelihood function of the Slots-plus-Averaging model
% ------------
% LLH=Slots_plus_Averaging(param, Data, Input)
%
% ## Theory ##
% This model assumed that there are a limited number of resouce slots and
% the fidelity of each representation depends on the number of slots being
% allocated to the corresponding item.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% - param
% K, kappa_1, kappa_r, (bias, biasF, s)
% - Data
% Data.error (response-sample), Data.SS (set size)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt)
% - Input
% Input.Variants (options of Variants), Input.PDF
%
% ## Output ##
% Summed log Likelihood
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
% Under the guidance of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function LLH=Slots_plus_Averaging(param, Data, Input)

% Specify parameters
K=param(1); % Capacity
kappa_1=param(2); % Unit precision
kappa_r=param(3); % Response variability
Nparam=3;
if Input.Variants.Bias==0
    bias=0; % Responses concentrate on samples
else
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
end
if Input.Variants.BiasF==0
    biasF=0; % Set bias as a consistent value
else
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
end
if Input.Variants.PrecF==0
    precF=0; % Set precision as consistent within each set size
else
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
end
if Input.Variants.Swap==0
    s=0; % No swap
else
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
end

% Configuration
errors=Data.error;
error_range=Data.error_range;
if length(error_range)==2
    continuous=1;
else
    continuous=0;
end
SS=Data.SS;
SS_range=unique(SS);
K=floor(K); % Note that capacity is a fixed, discrete value here
if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
    sample_range=Data.sample_range;
    samples=Data.sample;
else
    samples=ones(1,length(errors));
    sample_range=1;
end
if Input.Variants.Swap==1
    errors_nt=Data.error_nt;
    if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
        samples_nt=Data.sample_nt;
    end
end

% LH function
if continuous==1
    p_error=zeros(1,length(errors));
    p_error_NT=zeros(1,length(errors));
    kappa0=exp(log(kappa_1)*(cosd(4*samples)).^precF); % Flucutuative precision
    for i_error=1:length(errors)
        N=SS(i_error);
        kappa=kappa0(i_error);
        error0=errors(i_error)+bias+biasF*cosd(4*samples(i_error)-90);
        p_low=1-mod(K,N)/N;
        p_high=mod(K,N)/N;
        kappa_high=kappa*(floor(K/N)+1);
        kappa_low=kappa*floor(K/N);
        if K<N % Beyond capacity
            conv_1=sqrt((kappa).^2+(kappa_r)^2+2*kappa*kappa_r.*cosd(error0));
            p_error(i_error)=(1-K/N)*1/(error_range(2)-error_range(1))+...
                (K/N)*(besseli0_fast(conv_1)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r)));
        else
            conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0));
            conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0));
            p_error(i_error)=p_low*(besseli0_fast(conv_low)./(2*pi*besseli0_fast(kappa_low)*besseli0_fast(kappa_r)))+...
                p_high*(besseli0_fast(conv_high)./(2*pi*besseli0_fast(kappa_high)*besseli0_fast(kappa_r)));
        end
        if Input.Variants.Swap==1
            if N==1
                p_error_NT(i_error)=0;
            else
                p_temp_NT=0;
                for i_nt=1:N-1
                    error0_nt=errors_nt(i_nt, i_error)+bias+biasF*cosd(4*samples(i_error)-90); % Errors with bias
                    if K<N % beyond capacity
                        conv_1=sqrt((kappa).^2+(kappa_r)^2+2*kappa*kappa_r.*cosd(error0_nt));
                        p_temp_NT=p_temp_NT+(1-K/N)*1/(error_range(2)-error_range(1))+...
                            (K/N)*(besseli0_fast(conv_1)./(2*pi*besseli0_fast(kappa)*besseli0_fast(kappa_r)));
                    else % within capacity
                        conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0_nt));
                        conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0_nt));
                        p_temp_NT=p_temp_NT+p_low*(besseli0_fast(conv_low)./(2*pi*besseli0_fast(kappa_low)*besseli0_fast(kappa_r)))+...
                            p_high*(besseli0_fast(conv_high)./(2*pi*besseli0_fast(kappa_high)*besseli0_fast(kappa_r)));
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
    bias_cur=bias+biasF*cosd(4*sample_range-90); % Current bias
    kappa0=exp(log(kappa_1)*(cosd(4*sample_range)).^precF); % Current precision
    
    p_error=zeros(length(SS_range), length(error_range),length(sample_range));
    for i_N=1:length(SS_range)
        N=SS_range(i_N);
        if K<N % Beyond capacity
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias_cur;
                conv_1=sqrt((kappa0).^2+(kappa_r)^2+2*kappa0*kappa_r.*cosd(error0));
                p_error(i_N,i_error,:)=(1-K/N)*1/length(error_range)+...
                    (K/N)*(besseli0_fast(conv_1)./(2*pi*besseli0_fast(kappa0)*besseli0_fast(kappa_r)));
            end
            for i_sample=1:length(sample_range)
                p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
            end
        else
            p_low=1-mod(K,N)/N;
            p_high=mod(K,N)/N;
            kappa_high=kappa0*(floor(K/N)+1);
            kappa_low=kappa0*floor(K/N);
            for i_error=1:length(error_range)
                error0=error_range(i_error)+bias_cur;
                conv_low=sqrt((kappa_low).^2+(kappa_r)^2+2*kappa_low*kappa_r.*cosd(error0));
                conv_high=sqrt((kappa_high).^2+(kappa_r)^2+2*kappa_high*kappa_r.*cosd(error0));
                p_error(i_N,i_error,:)=p_low*(besseli0_fast(conv_low)./(2*pi*besseli0_fast(kappa_low)*besseli0_fast(kappa_r)))+...
                    p_high*(besseli0_fast(conv_high)./(2*pi*besseli0_fast(kappa_high)*besseli0_fast(kappa_r)));
            end
            for i_sample=1:length(sample_range)
                p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
            end
        end
    end
    
    % Calculate LH
    p_T=zeros(1,length(errors));
    p_NT=zeros(1,length(errors));
    for i=1:length(errors)
        if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
            p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), sample_range==samples(i));
        else
            p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), 1);
        end
    end
    if Input.Variants.Swap==1
        for i=1:length(errors_nt)
            if SS(i)==1
                p_NT(i)=0;
            else
                for j=1:SS(i)-1
                    if Input.Variants.BiasF==1 || Input.Variants.PrecF==1
                        p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(j,i), sample_range==samples_nt(j,i));
                    else
                        p_NT(i)=p_NT(i)+s/(SS(i)-1)*p_error(SS_range==SS(i),error_range==errors_nt(j,i), 1);
                    end
                end
            end
        end
    end
    p_LH=p_T+p_NT; % Target + non-target LH
end

% LLH
if isfield(Input,'PDF') && Input.PDF==1
    if Input.Variants.PrecF==1 || Input.Variants.BiasF==1
        if length(unique(Data.sample_range))~=1
            error('PDF mode requires sample to be unique')
        end
    elseif length(unique(Data.SS))~=1
        error('PDF mode requires set size to be unique')
    else
        LLH=p_error; % PDF
    end
else
    LLH=-sum(log(p_LH)); % Negative LLH
    if abs(LLH)==Inf || isnan(LLH)
        LLH=exp(666); % LLH should be a real value
    end
end

end
