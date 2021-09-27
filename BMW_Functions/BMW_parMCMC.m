%% Monte Carlo Markov Chain
%
% Sampling by Monte Carlo Markov Chain Algorithm (Adaptive Metropolis-Hastings/Differential Evolution)
%
% -----------------------
% [RawSampling, Summary]=BMW_parMCMC(model, data, config)
%
% ## Input ##
% - model
%   struct, contains the name of the model of interest
%   Type BMW('models') for the list of available models.
%       ~.Model
%           string, contains the name of the model of interest
%       ~.Constraints
%           array, contains ~.start(start values), ~.ub(upper bounds), ~.lb(lower bounds).
%           Note that the number of rows of model.Constraints.start stands for
%           the number of markov chains, each chain will non-repetitively
%           pick a row as its start value
%       ~.Variants
%           e.g. model.Variants.Bias==1;
%           See manual for details.
% - data
%   mat, the data variable that qualifies the simulation.
%   Type BMW('manual') for detailed requirements.
%
% - config
%   config.Algorithm
%       string, 'DE'/'MH', use 'DE' as default
%   config.MCMCparam
%       vector, configures the sampling
%       If we use 'DE', then the length of the vector should be two
%       and the elements stand for [gamma, eps] in order.
%           gamma, scaling factor of the between-chain difference vector
%           eps, the range of the uniform jitter
%       If we use 'MH', then the length should be three and the elements
%       should stand for [Sd, t0] in order.
%   config.Convergence
%       ~.Diagnostic
%           string, controls the method of convergence diagnostic, default as 'GR'
%           Available methods include 'GR' & 'wcGR'
%       ~.Nbatchburnin
%           integer, designates the number of burn-in samples in each batch
%           of sampling
%       ~.Nmaxbatchburnin
%           integer, designates the max quantity of burn-in batches
%   config.Verbosity
%       string, 'off'/'iter'/'notify'/'final', default as 'iter'
%       controls the level of display during sampling
%   config.GetParam
%       string, 'mean'/'mean_trimmed'/'median'/'median_trimmed'/'max'
%       denotes the method of finding the best parameter vector(sample)
%
% ## Output ##
% - RawSampling
% - Summary
%
% ## Reference ##
% -
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 3/3/2020
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [RawSampling, Summary]=BMW_parMCMC(model, data, config)

%% Prologue
% set default
Nparam=size(model.Constraints.start,2); % number of dimensions
Nchain=size(model.Constraints.start,1); % number of chains
if ~isfield(config,'Algorithm')
    config.Algorithm='DE'; % set DE as default algorithm 
end
if ~isfield(config,'Nsample')
    config.Nsample=max(5000, 2000*Nparam); % set default number of samples after convergence
end
if ~isfield(config,'Verbosity')
    config.Verbosity='iter';
end
Verbosity=config.Verbosity;
if ~isfield(config,'Transform')
    config.Transform='Probit';
end
if ~isfield(config,'Convergence')
    config.Convergence.Diagnostic='GR'; % set conventional GR as the default way to diagnose convergence
    config.Convergence.Nbatchburnin=200; % number of burn-in samples per batch
    config.Convergence.Nmaxbatchburnin=50; % max number of burn-in batches
    config.Convergence.Tol=0.1; % R threshold
else
    if ~isfield(config.Convergence,'Diagnostic'), config.Convergence.Diagnostic='GR'; end
    if ~isfield(config.Convergence,'Nbatchburnin'), config.Convergence.Nbatchburnin=200; end
    if ~isfield(config.Convergence,'Nmaxbatchburnin'), config.Convergence.Nmaxbatchburnin=50; end
    if ~isfield(config.Convergence,'Tol'), config.Convergence.Tol=0.1; end%+log(Nparam)/10; end
end
if ~isfield(config,'Ncore')
    config.Ncore=Nchain;
end
if ~isfield(config,'AutoParallel')
    config.AutoParallel=1;
end
if ~isfield(config,'GetParam')
    config.GetParam='mean_trimmed';
end
if strcmp(config.Transform,'NoTransform')
    Transform=0;
else
    Transform=1;
end
% set initial values
Nburnin=config.Convergence.Nbatchburnin; % number of burn-in samples per batch
MAXbatchburnin=config.Convergence.Nmaxbatchburnin; % max number of burn-in batches per chain
Nstate=ceil(config.Nsample/Nchain); % number of samples after convergence per chain
Nbatchburnin=Nburnin; % number of burn-in samples per batch per chain
model.RealOutput=model.Output; % record real output mode (LLH/LP/Prior/LPPD)
model.Output='All';
% set tuning parameters
if strcmp(config.Algorithm,'DE') % Differential Evolution
    if ~isfield(config,'MCMCparam') % default
        % Note that here we adopted the proposed values in Turner et al., 2013,
        % which might partly depend on Gelman et al., 1996,
        % with some adaptation
        gamma=2.38.^2/Nparam; % scales the difference vector
        eps=0.001; % random jitter (the range of the uniform distribution)
    else
        gamma=config.MCMCparam(1);
        eps=config.MCMCparam(2:end);
    end
elseif strcmp(config.Algorithm,'MH') % Metropolis-Hastings
    % We adopt the proposed defaults in Haario, Saksman, & Tamminen, 2001
    % which might depend on Gelman et al., 1996 as well.
    if Transform==0
        cov=((min(model.Constraints.ub-model.Constraints.lb)/10)^2/Nparam)*eye(Nparam); % initial covariance matrix of the multivariate gaussian transition kernel
    else
        cov=(2.38^2/Nparam)*eye(Nparam);
    end
    if ~isfield(config,'MCMCparam') % default
        Sd=2.38^2/Nparam; % scales the estimated covariance of the past samples
        t0=2*Nparam; % the start point of adaptation, assign 0 to call off the adaptation
        eps=0.001; % scales the constant that avoids the covariance matrix from being singular
    else
        Sd=config.MCMCparam(1);
        t0=config.MCMCparam(2);
        eps=config.MCMCparam(3);
    end
else
    warning('The algorithm is invalid, use DE-MCMC instead.')
    config.Algorithm='DE';
end
% get Ntrial
NewFunction=str2func(model.Model);
TestLH=NewFunction(model.Constraints.start, data, model);
Ntrial=length(TestLH.LPPD);
% pre-allocation
TrueSamples=zeros(Nchain, Nparam, Nstate);
Samples=zeros(Nchain,Nparam,Nstate);
logPosterior=zeros(Nchain,Nstate);
logPrior=zeros(Nchain,Nstate);
logLikelihood=zeros(Nchain,Nstate);
logPointwiseLH=zeros(Nchain,Ntrial,Nstate);
ConvergenceStat=zeros(1,MAXbatchburnin);
% shuffle random seed
rng(now);
% check parpool
Ncore=config.Ncore;
if Ncore>1 && rem(Ncore,1)==0
    if isempty(gcp('nocreate'))
        if config.AutoParallel==0
            MoveOn=input('\nParellel computing has not been activated yet, do you want to activate it now? (y/n) \n','s');
            fprintf('\n');
        elseif config.AutoParallel==1
            MoveOn='y';
        end
        if strcmp(MoveOn,'y')
            parpool([1,Ncore],'AttachedFiles',{model.Model});
        else
            fprintf('\Nstate right. Run MCMC with the single worker...\n')
        end
    end
elseif Ncore<1 || rem(Ncore,1)~=0
    error('The number of cores is invalid...')
end

%% MCMC
if strcmp(config.Algorithm,'DE') % Differential Evolution
    % get start value
    start=MCMCConvert_BMW(model.Constraints.start,model.Constraints.ub,model.Constraints.lb,config.Transform);
    % burn-in sampling
    fprintf('\nNow start DE-MCMC sampling...\n')
    for t_burnin=1:MAXbatchburnin
        BS=DEMCMCchain([gamma, eps],Nbatchburnin,Nchain,config.Transform,start,data,model);
        start=BS.Samples(:,:,end);
        BS.TrueSamples(:,:,end)
        [BoolConvergence, ConvergenceStat(t_burnin)]=TestConvergence(BS.TrueSamples,config.Convergence.Tol,config.Convergence.Diagnostic);
        if BoolConvergence
            if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
                fprintf('\n%d burn-in samples collected. R=%d\n',t_burnin*Nbatchburnin*Nchain,ConvergenceStat(t_burnin))
                fprintf('\nConverged! Now start collecting valid samples...\n')
            end
            break;
        else
            if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
                fprintf('\n%d burn-in samples collected. R=%d\n',t_burnin*Nbatchburnin*Nchain,ConvergenceStat(t_burnin))
            end
        end
    end
    % collect valid samples
    fprintf('\n')
    Nstate_rest=Nstate;
    for batch=1:10
        Nstate_cur=min(ceil(Nstate/10),Nstate_rest);
        S=DEMCMCchain([gamma, eps],Nstate_cur,Nchain,config.Transform,start,data,model);
        start=S.Samples(:,:,end);
        Samples(:,:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.Samples;
        TrueSamples(:,:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.TrueSamples;
        logPosterior(:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.LP;
        logLikelihood(:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.LLH;
        logPrior(:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.Prior;
        logPointwiseLH(:,:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.LPPD;
        Nstate_rest=Nstate_rest-Nstate_cur;
        if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
            fprintf('|||||')
        end
    end
    if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
        fprintf('\n\nDone! %d samples collected in total.\n',Nstate*Nchain)
    end
    
elseif strcmp(config.Algorithm,'MH') % Metropolis-Hastings
    % get start value
    start=MCMCConvert_BMW(model.Constraints.start,model.Constraints.ub,model.Constraints.lb,config.Transform);
    % burn-in sampling
    for t_burnin=1:MAXbatchburnin
        [BS,cov]=MHMCMCchain([Sd, t0, eps],cov,Nbatchburnin,Nchain,config.Transform,start,model);
        start=BS.Samples(:,:,end);
        [BoolConvergence, ConvergenceStat(t_burnin)]=TestConvergence(BS.TrueSamples,config.Convergence.Tol,config.Convergence.Diagnostic);
        if BoolConvergence
            if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
                fprintf('\nConverged! Now start collecting valid samples...\n')
            end
            break;
        else
            if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
                fprintf('\n%d burn-in samples collected. R=%d\n',t_burnin*Nbatchburnin*Nchain,ConvergenceStat(t_burnin))
            end
        end
    end
    % collect valid samples
    fprintf('\n')
    Nstate_rest=Nstate;
    for batch=1:10
        Nstate_cur=min(ceil(Nstate/10),Nstate_rest);
        [S,cov]=MHMCMCchain([Sd, t0, eps],cov,Nstate_cur,Nchain,config.Transform,start,data,model);
        start=S.Samples(:,:,end);
        Samples(:,:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.Samples;
        TrueSamples(:,:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.TrueSamples;
        logPosterior(:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.LP;
        logLikelihood(:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.LLH;
        logPrior(:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.Prior;
        logPointwiseLH(:,:,Nstate-Nstate_rest+1:Nstate-Nstate_rest+Nstate_cur)=S.LPPD;
        Nstate_rest=Nstate_rest-Nstate_cur;
        if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
            fprintf('|||||')
        end
    end
    if strcmp(Verbosity,'iter') || strcmp(Verbosity,'detail')
        fprintf('\n\nDone! %d samples collected in total.\n',Nstate*Nchain)
    end
    
end

%% Merge Chains
if Transform==0, TrueSamples=Samples; end
RawSampling.RawSamples=zeros(Nstate*Nchain, Nparam);
RawSampling.Samples=zeros(Nstate*Nchain, Nparam);
RawSampling.logPosterior=zeros(Nstate*Nchain, 1);
RawSampling.logLikelihood=zeros(Nstate*Nchain, 1);
RawSampling.logPointwiseLikelihood=zeros(Nstate*Nchain, Ntrial);
RawSampling.logPrior=zeros(Nstate*Nchain, 1);
RawSampling.Chain=zeros(Nstate*Nchain, 1);
for chain=1:Nchain
    if Transform==1
        RawSampling.RawSamples((chain-1)*Nstate+1:chain*Nstate,:)=permute(Samples(chain, :, :),[3, 2, 1]); % raw sample values
    end
    RawSampling.Samples((chain-1)*Nstate+1:chain*Nstate,:)=permute(TrueSamples(chain, :, :),[3, 2, 1]); % sample values
    RawSampling.logPosterior((chain-1)*Nstate+1:chain*Nstate)=logPosterior(chain, :)'; % Posterior values
    RawSampling.logLikelihood((chain-1)*Nstate+1:chain*Nstate)=logLikelihood(chain, :)'; % Posterior values
    RawSampling.logPrior((chain-1)*Nstate+1:chain*Nstate)=logPrior(chain, :)'; % Posterior values
    RawSampling.logPointwiseLikelihood((chain-1)*Nstate+1:chain*Nstate,:)=permute(logPointwiseLH(chain, :, :),[3, 2, 1]); % Posterior values
    RawSampling.Chain((chain-1)*Nstate+1:chain*Nstate)=chain*ones(Nstate,1); % chain indices
end

%% Summary Statistics
model.Output=model.RealOutput;
switch config.GetParam
    case 'max'
        [Summary.BestPosterior, IndMAXposterior]=max(RawSampling.logPosterior); % max of logPosterior density
        Summary.FitParam=RawSampling.Samples(IndMAXposterior,:); % best parameter(s)
    case 'mean'
        Summary.FitParam=mean(RawSampling.Samples); % best parameter(s)
        Summary.BestPosterior=NewFunction(Summary.FitParam, data, model);
    case 'mean_trimmed'
        Samples_Trimmed=sort([RawSampling.logPosterior,RawSampling.Samples]);
        min_Trim=ceil(0.05/2*size(Samples_Trimmed,1));
        max_Trim=size(Samples_Trimmed,1)-ceil(0.05/2*size(Samples_Trimmed,1));
        Summary.FitParam=mean(Samples_Trimmed(min_Trim:max_Trim,2:end)); % best parameter(s)
        Summary.BestPosterior=NewFunction(Summary.FitParam, data, model);
    case 'median'
        Summary.FitParam=median(RawSampling.Samples); % best parameter(s)
        Summary.BestPosterior=NewFunction(Summary.FitParam, data, model);
    case 'median_trimmed'
        Samples_Trimmed=sort([RawSampling.logPosterior,RawSampling.Samples]);
        min_Trim=ceil(0.05/2*size(Samples_Trimmed,1));
        max_Trim=size(Samples_Trimmed,1)-ceil(0.05/2*size(Samples_Trimmed,1));
        Summary.FitParam=median(Samples_Trimmed(min_Trim:max_Trim,2:end)); % best parameter(s)
        Summary.BestPosterior=NewFunction(Summary.FitParam, data, model);
end

%% Epilogue
if ~strcmp(config.Verbosity,'off')
    fprintf('All done!\n\n')
end

end

%% DE-MCMC Chains
function RawSamples=DEMCMCchain(Param,Nstate,Nchain,Transform,start,data,model)
% define parameters
gamma=Param(1);
eps=Param(2);
Constraints=model.Constraints;
NewFunction=str2func(model.Model);
TestLH=NewFunction(MCMCConvert_BMW(start(1,:),model.Constraints.ub,model.Constraints.lb,['Inverse', Transform]), data, model);
Ntrial=length(TestLH.LPPD);
start1=start;
% run chains
LP=zeros(Nchain,Nstate);
LLH=zeros(Nchain,Nstate);
Prior=zeros(Nchain,Nstate);
LPPD=zeros(Nchain,Ntrial,Nstate);
Samples=zeros(Nchain,size(start,2),Nstate);
TrueSamples=zeros(Nchain,size(start,2),Nstate);
for state=1:Nstate
    parfor chain=1:Nchain
        [LP(chain,state),LLH(chain,state),Prior(chain,state),LPPD(chain,:,state),Samples(chain,:,state),TrueSamples(chain,:,state)]=...
            DEMCMCchain_single(Transform,chain,Nchain,start1,NewFunction,gamma,eps,Constraints,model,data);
    end
    start1=Samples(:,:,state);
end
RawSamples.LP=LP;
RawSamples.LLH=LLH;
RawSamples.Prior=Prior;
RawSamples.LPPD=LPPD;
RawSamples.Samples=Samples;
RawSamples.TrueSamples=TrueSamples;

end

function [LP,LLH,Prior,LPPD,Samples,TrueSamples]=DEMCMCchain_single(Transform,chain,Nchain,start1,NewFunction,gamma,eps,Constraints,model,data)
% generate proposal state
ChainRange=1:Nchain; % decision space of chains
ChainRange=ChainRange(ChainRange~=chain); % exclude the current chain
while 1 % redo sampling if the current sample is out of bounds
    ChainLot=randsample(ChainRange,2); % ramdomly pick two chains without repetition
    ChainDiff=gamma*(start1(ChainLot(1),:)-start1(ChainLot(2),:)); % get difference vector
    DiffNum=randsample(1:size(start1,2),1);
    ChainDiff(DiffNum)=randsample([0, 5, ones(1,8)],1)*ChainDiff(DiffNum);
    UniNoise=2*eps.*rand(1,length(eps))-abs(eps); % noise
    PropState=start1(chain,:)+ChainDiff+UniNoise; % proposal state
    if strcmp(Transform,'NoTransform')
        if ~any(PropState>Constraints.ub) && ~any(PropState<Constraints.lb) % test bounds
            break;
        end
    else
        break;
    end
end
% test proposal state
TruePropState=MCMCConvert_BMW(PropState,Constraints.ub,Constraints.lb,['Inverse', Transform]);
TruePrevState=MCMCConvert_BMW(start1(chain,:),Constraints.ub,Constraints.lb,['Inverse', Transform]);
minMHr=rand; % sample Metropolis-Hastings ratio threshold
Prop0=NewFunction(TruePropState, data, model);
Prev0=NewFunction(TruePrevState, data, model);
%             eval(['Prop0=',Model,'(TruePropState, data, model);'])
%             eval(['Prev0=',Model,'(TruePrevState, data, model);'])
PropPosterior=-Prop0.(model.RealOutput);
PrevPosterior=-Prev0.(model.RealOutput);
%             eval(['PropPosterior=-Prop0.',RealOutput,';'])
%             eval(['PrevPosterior=-Prev0.',RealOutput,';'])
MHr=exp(PropPosterior-PrevPosterior); % current MH ratio
if MHr>minMHr
    CurrentState=PropState; % accept proposal
    if strcmp(Transform,'NoTransform')
        TrueCurrentState=PropState;
    else
        TrueCurrentState=TruePropState;
    end
    LP=PropPosterior; % record log posterior
    LLH=-Prop0.LLH; % record log likelihood
    Prior=-Prop0.Prior; % record prior
    LPPD=-Prop0.LPPD; % record pointwise likelihood
else
    CurrentState=start1(chain,:); % reject, use the previous state instead
    if strcmp(Transform,'NoTransform')
        TrueCurrentState=start1(chain,:);
    else
        TrueCurrentState=TruePrevState;
    end
    LP=PrevPosterior; % record log posterior
    LLH=-Prev0.LLH; % record log likelihood
    Prior=-log(Prev0.Prior); % record prior
    LPPD=-Prev0.LPPD; % record pointwise likelihood
end
Samples=CurrentState; % record sample value
TrueSamples=TrueCurrentState;
end

%% MH-MCMC Chains
function [RawSamples,cov]=MHMCMCchain(Param,cov0,Nstate,Nchain,Transform,start,data,model)

% define parameters
Sd=Param(1);
t0=Param(2);
eps=Param(3);
Constraints=model.Constraints;
NewFunction=str2func(model.Model);
TestLH=NewFunction(MCMCConvert_BMW(start(1,:),model.Constraints.ub,model.Constraints.lb,['Inverse', Transform]), data, model);
Ntrial=length(TestLH.LPPD);
start1=start;
% run chains
AllOutputs=cel(1,Nchain);
for i=1:Nchain
    AllOutputs{i}.LP=zeros(1,Nstate);
    AllOutputs{i}.LLH=zeros(1,Nstate);
    AllOutputs{i}.Prior=zeros(1,Nstate);
    AllOutputs{i}.LPPD=zeros(1,Ntrial,Nstate);
    AllOutputs{i}.Samples=zeros(1,size(start,2),Nstate);
    AllOutputs{i}.TrueSamples=zeros(1,size(start,2),Nstate);
end
parfor chain=1:Nchain
    [AllOutputs{chain},cov{chain}]=...
        MHMCMCchain_single(Transform,Nstate,start1(chain,:),NewFunction,cov0,Sd,t0,eps,Constraints,model,data);
    %[LP(chain),LLH(chain,state),Prior(chain,state),LPPD(chain,:,state),Samples(chain,:,state),TrueSamples(chain,:,state)]=...
    %  DEMCMCchain_single(Transform,chain,Nchain,start1,NewFunction,Sd,t0,eps,Constraints,model,data);
end

% Combine chains
LP=zeros(Nchain,Nstate);
LLH=zeros(Nchain,Nstate);
Prior=zeros(Nchain,Nstate);
LPPD=zeros(Nchain,Ntrial,Nstate);
Samples=zeros(Nchain,size(start,2),Nstate);
TrueSamples=zeros(Nchain,size(start,2),Nstate);
for i=1:Nchain
    LP(i,:)=AllOutputs{i}.LP;
    LLH(i,:)=AllOutputs{i}.LLH;
    Prior(i,:)=AllOutputs{i}.Prior;
    LPPD(i,:,:)=AllOutputs{i}.LPPD;
    Samples(i,:,:)=AllOutputs{i}.Samples;
    TrueSamples(i,:,:)=AllOutputs{i}.TrueSamples;
end

RawSamples.LP=LP;
RawSamples.LLH=LLH;
RawSamples.Prior=Prior;
RawSamples.LPPD=LPPD;
RawSamples.Samples=Samples;
RawSamples.TrueSamples=TrueSamples;

end

function [AllOutputs,cov]=MHMCMCchain_single(Transform,Nstate,start1,NewFunction,cov,Sd,t0,eps,Constraints,model,data)

for state=1:Nstate
    % update proposal distribution
    if state==1
        History=start1; % start
    else
        % append the past samples to the start value
        % use inverse to avoid the problem of Nparam and (state-1) being equal
        History=vertcat(start1,reshape(AllOutputs.Samples(1,:,1:(state-1)),[Nparam, state-1])');
    end
    if mod(state,10)==0
        cov=UpdateCov(cov,History,Sd,t0,eps); % adapt the tuning parameters based on the history of sampling
    end
    while 1
        % generate proposal state
        PrevState=History(end,:); % get previous state
        PropState=mvnrnd(PrevState,cov); % sample from multivariate normal distribution
        if strcmp(Transform,'NoTransform')
            if ~any(PropState>Constraints.ub) && ~any(PropState<Constraints.lb) % test bounds
                break;
            end
        else
            break;
        end
    end
    TruePropState=MCMCConvert_BMW(PropState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
    TruePrevState=MCMCConvert_BMW(PrevState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
    % test proposal state
    minMHr=rand; % sample Metropolis-Hastings ratio threshold
    Prop0=NewFunction(TruePropState, data, model);
    Prev0=NewFunction(TruePrevState, data, model);
    %eval(['Prop0=',model.Model,'(TruePropState, data, model);'])
    %eval(['Prev0=',model.Model,'(TruePrevState, data, model);'])
    PropPosterior=-Prop0.(model.RealOutput);
    PrevPosterior=-Prev0.(model.RealOutput);
    %eval(['PropPosterior=-Prop0.',RealOutput,';'])
    %eval(['PrevPosterior=-Prev0.',RealOutput,';'])
    MHr=exp(PropPosterior-PrevPosterior)*(mvnpdf(PrevState,PropState,cov)/mvnpdf(PropState,PrevState,cov)); % current MH ratio
    if MHr>minMHr
        CurrentState=PropState; % accept proposal
        if strcmp(Transform,'NoTransform')
            TrueCurrentState=PropState;
        else
            TrueCurrentState=TruePropState;
        end
        AllOutputs.LP(1,state)=PropPosterior; % record logPosterior
        AllOutputs.LLH(1,state)=-Prop0.LLH; % record likelihood
        AllOutputs.Prior(1,state)=-Prop0.Prior; % record prior;
        AllOutputs.LPPD(1,:,state)=-Prop0.LPPD; % record pointwise likelihood
    else
        CurrentState=PrevState; % reject, use the previous state instead
        if strcmp(Transform,'NoTransform')
            TrueCurrentState=PrevState;
        else
            TrueCurrentState=TruePrevState;
        end
        AllOutputs.LP(1,state)=PrevPosterior; % record logPosterior
        AllOutputs.LLH(1,state)=-Prev0.LLH; % record likelihood
        AllOutputs.Prior(1,state)=-log(Prev0.Prior); % record prior;
        AllOutputs.LPPD(1,:,state)=-Prev0.LPPD; % record pointwise likelihood
    end
    AllOutputs.Samples(1,:,state)=CurrentState; % record sample value
    if Transform==1, AllOutputs.TrueSamples(1,:,state)=TrueCurrentState; end
end

end

%% Update the Tuning Parameter of Adaptive MH-MCMC
% Based on the adaptive MCMC algorithm in Haario, Saksman, & Tamminen, 2001
function C1=UpdateCov(C0, History, Sd, t0, eps)
if t0==0, t0=Inf; end % no adaptation
Nparam=size(History,1); % number of dimensions
t=size(History,2); % number of past steps
if t<=t0 % start
    C1=C0; % use the initial parameters
else % update cov matrix based on the sampling history
    % Note that the history here contains the start values
    MeanSample0=sum(History(1:t-1),2)/(t-1);
    MeanSample=sum(History,2)/t;
    CurrentSample=History(:,end);
    % get empirical estimated covariance of the past samples
    C1=(t-2)*C0/(t-1)+(Sd/(t-1))*((t-1)*(MeanSample0)'*MeanSample0-...
        t*(MeanSample)'*MeanSample+CurrentSample'*CurrentSample+...
        eps*eye(Nparam)); % add this constant to ensure boundary convergence
end
end


%% Diagnose Convergence for Iterative Simulation
% Generally refer to Roy, 2019
function [BoolConverge, Stat]=TestConvergence(History, config, method)
if strcmp(method,'GR'), [BoolConverge, Stat]=TestConvergence_GR(History, config);...
elseif strcmp(method,'wcGR'), [BoolConverge, Stat]=TestConvergence_wcGR(History, config);...
end
end

% TestConvergence_GR based on GR statistic in Gelman & Rubin, 1992
function [BoolConverge, R]=TestConvergence_GR(History, config)
% set default
delta=config;
Nsample=size(History, 3); % number of samples
Nchain=size(History, 1); % number of chains
Smean=mean(History, 3); % mean across sample (Nchain*Nparam mat)
Allmean=mean(Smean, 1); % grand mean (Nparam vector)
Svar=sum((History-repmat(Smean,[1,1,Nsample])).^2, 3)/(Nsample-1); % mean sample variance
Svarchain=mean(Svar, 1); % mean sample variance across chain
B=sum((Smean-repmat(Allmean,[Nchain,1])).^2)/(Nchain-1); % between-chain variance
EstVar=Svarchain*(Nsample-1)/Nsample+B; % estimated variance of target distribution
R=sqrt(EstVar/Svarchain); % GR potential scale reduction factor (PSRF)
if R<1+delta
    BoolConverge=1; % converged
else
    BoolConverge=0; % not coverged yet
end
end

% TestConvergence_wcGR based on the adapted GR statistic in Vats & Knudson, 2019
function [BoolConverge, R]=TestConvergence_wcGR(History, config)
% set default
delta=config;
Nsample=size(History, 3);
Smean=mean(History, 3); % mean across sample (Nchain*Nparam mat)
% decide batch size, we use sqrt(Nsample) here
B1=floor(sqrt(Nsample)); A1=floor(Nsample/B1);
B2=B1/3; A2=floor(Nsample/B1);
% define the lugsail batch means estimator
    function tau=lugsail(a,b,X,m)
        Y=zeros(size(X,1),size(X,2),a); % Nchain*Nparam*Nbatch mat
        for batch=1:a % loop batch
            Y(:,:,batch)=sum(X(:,:,(batch-1)*b+1:batch*b),3)/b;
        end
        tau=sum((Y-repmat(m,[1,1,a])).^2, 3)*b/(a-1);
    end
tau1=lugsail(A1,B1,History,Smean);
tau2=lugsail(A2,B2,History,Smean);
tau_all=mean(2*tau1-tau2,1); % the final estimator
Svar=sum((History-repmat(Smean,[1,1,Nsample])).^2, 3)/(Nsample-1); % mean sample variance
Svarchain=mean(Svar, 1); % mean sample variance across chain
EstVar=Svarchain*(Nsample-1)/Nsample+tau_all/Nsample; % estimated variance of target distribution
R=sqrt(EstVar/Svarchain); % GR potential scale reduction factor (PSRF)
if R<1+delta
    BoolConverge=1; % converged
else
    BoolConverge=0; % not coverged yet
end

end
