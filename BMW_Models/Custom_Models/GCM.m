%% Generalized Context Model
%
% LLH function of Generalized Context Model(GCM)
% ------------
% ## Input ##
% - param
% w(:) (dimension weights), c (sensory scaling)
% (y (response scaling))
% - Data
% Data.exemplar (#exemplar-by-dimension), Data.exemplar_category (#category of exemplars)
% Data.response (reported category in each trial (#trial-by-dimension)),
% Data.sample (sample in each trial (#trial-by-dimension))
% - Input
% Input.Derivatives.Response (Response scaling)
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 11/5/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function LLH=GCM(param, Data, Input)
% Specify parameters
Nw=size(Data.sample,2); % Number of dimensions
w=param(1:Nw); % Weights
c=param(Nw+1); % Sensory scaling parameter
Nparam=Nw+1;
if isfield(Input,'Derivatives')
    if isfield(Input.Derivatives,'Response') && Input.Derivatives.Response==1;
        y=param(Nparam+1); % Response scaling parameter
    else
        y=1;
    end
end

% Configuration
samples=Data.sample;
sample_range=unique(samples,'rows');
E=Data.exemplar;
E_category=Data.exemplar_category; % Categories (1~N)
responses=Data.response; % 1~N
Nc=length(unique(E_category)); % # of categories

% Convert 'samples' into numbers
N_samples=zeros(length(samples),1);
for i=1:length(samples)
    for j=1:size(sample_range,1)
        N_samples(i)=N_samples(i)+j*isequal(samples(i,:),sample_range(j,:));
    end
end
% Distance & similarity
d=zeros(Nc,length(sample_range));
for cat=1:Nc
    E_temp=E(E_category==cat,:);
    for s=1:size(sample_range,1)
        for dim=1:Nw
            d(cat,s)=d(cat,s)+sum(w(dim)*abs(E_temp(:,dim)-sample_range(s,dim)*ones(size(E_temp,1),1)));
        end
    end
end
% Likelihood function
p=zeros(Nc,length(samples));
for s=1:size(sample_range,1)
    for cat=1:Nc
        p(cat,s)=(exp(-c*d(cat,s))).^y/sum((exp(-c*d(:,s))).^y);
    end
end

p_LH=zeros(1,length(responses));
for i=1:length(responses)
    p_LH(i)=p(responses(i),N_samples(i));
end

% LLH
LLH=-sum(log(p_LH)); % Negative LLH

if LLH==Inf || isnan(LLH)
    LLH=exp(666); % LLH should be a real value
end

end