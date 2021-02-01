%% Do simulation for resource-based WM theoretical models
%
% Simulate resource-based WM models
% config
%   .

function data_post=SimulationWM_BMW(param, data_prior, config)
data_post=data_prior;
SS_set=sort(unique(data_prior.SS));
Nss=length(SS_set);
Ntrial=size(data_prior.sample,1);
Nnt=size(data_prior.sample,2);
if nargin==2 || ~isfield(config,'model')
    config.model='Slots-plus-Averaging';
end
if nargin==2 || ~isfield(config,'BI')
    config.BI=[zeros(1,Nss), 1];
end
if nargin==2 || ~isfield(config,'WI')
    config.WI=[zeros(1,Nss), 1];
end
if nargin==2 || ~isfield(config,'NT')
    config.NT=0;
end
if nargin==2 || ~isfield(config,'Bias')
    config.Bias=0;
end
% discrete error or continuous error?
if length(data_prior.response_range)==2
    discrete=0;
else
    discrete=1;
end

Nmem=(Nnt-1)*double(config.NT~=0)+1;
% %% Capacity?
% if any(strcmp(config.model,{'Slots-plus-Averaging'}))
%     data_g=Sim_G(g, data_prior, Nmem, config);
% end
%% Continuous WM
switch config.model
    case 'Slots-plus-Averaging'
        [trial_type, data_o]=Sim_SA(param, data_prior, Nmem, config);
    case 'Variable Precision'
        [trial_type, data_o]=Sim_VP(param, data_prior, Nmem, config);
end

%% Between-Item Categorical WM
data_c=Sim_BI(data_prior, Nmem, config);

%% Combine
response=zeros(Ntrial,Nmem);
for n=1:Nmem
    for ss=1:Nss
        ss_trial_id=find(data_prior.SS==SS_set(ss));
        ng_trial_id=find(trial_type(n,:)>0);
        g_trial_id=find(trial_type(n,:)==0);
        ng_trial_id_all=sort(union(ng_trial_id, ss_trial_id));
        g_trial_id_all=sort(union(g_trial_id, ss_trial_id));
        Nng=length(ng_trial_id_all);
        Nctrial=binornd(Nng, config.BI(n));
        ctrial_id=sort(randsample(1:Nng,Nctrial));
        otrial_id=setdiff(1:Nng, ctrial_id);
        response_g=data_o.response(g_trial_id_all, n);
        response_o=data_o.response(ng_trial_id_all, n);
        response_c=data_c.response(ng_trial_id_all, n);
        response(g_trial_id_all,n)=response_g;
        response_ng=zeros(Nng,1);
        response_ng(ctrial_id)=response_c(ctrial_id);
        response_ng(otrial_id)=response_o(otrial_id);
        response(ng_trial_id_all,n)=response_ng;
    end
end
% discretize?
if discrete==1
    response_diff0=abs(CircDist_BMW('Diff',repmat(data_prior.response_range,[Ntrial,1]),repmat(response,[1,length(data_prior.response_range)])));
    [~,r_ind0]=min(response_diff0,[],2);
    response=data_prior.response_range(r_ind0)';
end
data_post.response=response;

end

%% Simulate SA
function [trial_type, data1]=Sim_SA(param, data0, Nmem, config)
K=param(1);
kappa_1=param(2);
SS_set=unique(data0.SS);
Nss=length(SS_set);
Ntrial=size(data0.sample,1);
Nnt=size(data0.sample,2);
w_c=config.WI(1:Nss);
kappa_c=config.WI(end);
response_range=data0.response_range;
% number of samples to generate / trial
if Nmem==1
    Nm=1;
else
    Nm=Nnt;
end
response=zeros(Ntrial,Nm);
trial_type=zeros(Ntrial,Nm);
% range
% if length(data0.error_range)==2
%     error_range=[data0.error_range(1), data0.error_range(2)];
% else
%     error_range=[min(data0.error_range), max(data0.error_range)];
% end
% simulation Lego!!!!
period=max(response_range)-min(response_range);
for t=1:Ntrial
    N=data0.SS(t);
    SS_id=find(SS_set==N);
    config_mem.weight=[1-w_c(SS_id), w_c(SS_id)];
    if N<K
        p_low=1-mod(K,N)/N;
        kappa_high=kappa_1*(floor(K/N)+1);
        kappa_low=kappa_1*floor(K/N);
        for n=1:Nm
            mu=CircSummary_BMW('Mean',[data0.sample(t,n),data0.category(t,n)],period,config_mem);
            mu=CircDist_BMW('Diff',mu,-config.Bias,period);
            if binornd(1,p_low)==1
                response_haha=vmrand(mu/180*pi,kappa_low,response_range/180*pi,1);
                trial_type(t,n)=1;
            else
                response_haha=vmrand(mu/180*pi,kappa_high,response_range/180*pi,1);
                trial_type(t,n)=2;
            end
            response_haha=response_haha/pi*180; % rad2ang
            response(t,n)=response_haha;
        end
    else
        for n=1:Nm
            if binornd(1,1-K/N)==1
                % generate guesses
                response(t,n)=(max(response_range)-min(response_range))*rand(1)-min(response_range);
                trial_type(t,n)=0;
            else
                % generate memory responses
                kappa=(1-w_c(SS_id))*kappa_1+w_c(SS_id)*kappa_c;
                mu=CircSummary_BMW('Mean',[data0.sample(t,n),data0.category(t,n)],period,config_mem);
                mu=CircDist_BMW('Diff',mu,-config.Bias,period);
                response_haha=vmrand(mu/180*pi,kappa,response_range/180*pi,1);
                response_haha=response_haha/pi*180; % rad2ang
                response(t,n)=response_haha;
                trial_type(t,n)=3;
            end
        end
    end
end
data1.response=response;
end

%% Simulate VP
function data1=Sim_VP(param, data0, Nmem, config)
K=param(1);
kappa_1=param(2);
SS_set=unique(data0.SS);
Nss=length(SS_set);
Ntrial=size(data0.sample,1);
Nnt=size(data0.sample,2);
w_c=config.WI(1:Nss);
kappa_c=config.WI(end);
%response_range=data0.response_range;
% number of samples to generate per trial
if Nmem==1
    Nm=1;
else
    Nm=Nnt;
end
response=zeros(Ntrial,Nm);
trial_type=zeros(Ntrial,Nm);
% range
if length(data0.error_range)==2
    error_range=[data0.error_range(1), data0.error_range(2)];
else
    error_range=[min(data0.error_range), max(data0.error_range)];
end
% simulation Lego!!!!
period=error_range(2)-error_range(1);
for t=1:Ntrial
    SS_id=find(SS_set==data0.SS(t));
    config_mem.weight=[1-w_c(SS_id), w_c(SS_id)];
    if data0.SS(t)<K
        p_low=1-mod(K,data0.SS(t))/data0.SS(t);
        kappa_high=kappa_1*(floor(K/data0.SS(t))+1);
        kappa_low=kappa_1*floor(K/data0.SS(t));
        for n=1:Nm
            mu=CircSummary_BMW('Mean',[data0.sample(t),data0.category(t,n)],period,config_mem);
            mu=CircDist_BMW('Diff',mu,-Bias,period);
            if binornd(1,p_low)==1
                response(t,n)=vmrand(mu,kappa_low,1);
                trial_type(t,n)=1;
            else
                response(t,n)=vmrand(mu,kappa_high,1);
                trial_type(t,n)=2;
            end
        end
    else
        for n=1:Nm
            if binornd(1,1-K/data0.SS(t))==1
                % generate guesses
                response(t,n)=(error_range(2)-error_range(1))*rand(1)-error_range(2);
                trial_type(t,n)=0;
            else
                % generate memory responses
                kappa=(1-w_c(SS_id))*kappa_1+w_c(SS_id)*kappa_c;
                mu=CircSummary_BMW('Mean',[data0.sample(t),data0.category(t,n)],period,config_mem);
                mu=CircDist_BMW('Diff',mu,-Bias,period);
                response(t,n)=vmrand(mu,kappa,1);
                trial_type(t,n)=3;
            end
        end
    end
end
data1.response=response;
end

%% Simulate BI
function data1=Sim_BI(data0, Nmem, config)
% SS_set=unique(data0.SS);
% Nss=length(SS_set);
Ntrial=size(data0.sample,1);
Nnt=size(data0.sample,2);
kappa_c=config.BI(end);
response_range=data0.response_range;
% number of samples to generate per trial
if Nmem==1
    Nm=1;
else
    Nm=Nnt;
end
response=zeros(Ntrial,Nm);
% range
% if length(data0.error_range)==2
%     error_range=[data0.error_range(1), data0.error_range(2)];
% else
%     error_range=[min(data0.error_range), max(data0.error_range)];
% end
% simulation Lego!!!!
period=max(response_range)-min(response_range);
% trial_type=zeros(Ntrial,Nm);
for t=1:Ntrial
    for n=1:Nm
        mu=data0.category(t,n);
        mu=CircDist_BMW('Diff',mu,-config.Bias,period);
        response_haha=vmrand(mu/180*pi,kappa_c,response_range/180*pi,1);
        response_haha=response_haha/pi*180; % rad2ang
        %response_haha=response_haha.*(max(response_range)-min(response_range))/360;
        response(t,n)=response_haha;
    end
end
data1.response=response;
end

%% Convert
% function data=convert(data, range)
% data=data/180*pi;
% while any(data<range(1)) || any(data>=range(2))
%     data(data>=range(2))=data(data>=range(2))-range(2);
%     data(data<range(1))=data(dasta<range(1))-range(1);
% end
% end
