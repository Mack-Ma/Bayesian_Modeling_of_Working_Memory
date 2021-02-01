%% Monte Carlo Markov Chain
%
% Sampling by Monte Carlo Markov Chain Algorithm (Adaptive Metropolis-Hastings/Differential Evolution)
%
% -----------------------
% [RawSampling, Summary]=BMW_MCMC(model, data, config)
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
%           Available methods include 'GR' & 'GRL'
%       ~.Nbatchburnin
%           integer, designates the number of burn-in samples in each batch
%           of sampling
%       ~.Nmaxbatchburnin
%           integer, designates the max quantity of burn-in batches
%   config.Verbosity
%       string, 'off'/'iter'/'notify'/'final', default as 'iter'
%       controls the level of display during sampling
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
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/BMW-Ma/Bayesian_Modeling_of_Working_Memory
%

function [RawSampling, Summary]=BMW_MCMC(model, data, config)

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
if ~isfield(config,'Transform')
    config.Transform='Probit';
end
if ~isfield(config,'Convergence')
    config.Convergence.Diagnostic='GR'; % set conventional GR as the default way to diagnose convergence
    config.Convergence.Nbatchburnin=200; % number of burn-in samples per batch
    config.Convergence.Nmaxbatchburnin=50; % max number of burn-in batches
    config.Convergence.Tol=0.1+log(Nparam)/10; % R threshold
else
    if ~isfield(config.Convergence,'Diagnostic'), config.Convergence.Diagnostic='GR'; end
    if ~isfield(config.Convergence,'Nbatchburnin'), config.Convergence.Nbatchburnin=200; end
    if ~isfield(config.Convergence,'Nmaxbatchburnin'), config.Convergence.Nmaxbatchburnin=25; end
    if ~isfield(config.Convergence,'Tol'), config.Convergence.Tol=0.1+log(Nparam)/10; end
end
if strcmp(config.Transform,'NoTransform')
    Transform=0;
else
    Transform=1;
end
% set initial values
Nburnin=config.Convergence.Nbatchburnin; % number of burn-in samples per batch per chain
MAXbatchburnin=config.Convergence.Nmaxbatchburnin; % max number of burn-in batches per chain
Nstate=ceil(config.Nsample/Nchain); % number of samples after convergence per chain
Nbatchburnin=ceil(Nburnin/Nchain); % number of burn-in samples per batch per chain
RealOutput=model.Output; % record real output mode (LLH/LP/Prior/LPPD)
model.Output='All';
% set tuning parameters
if strcmp(config.Algorithm,'DE') % Differential Evolution
    if ~isfield(config,'MCMCparam') % default
        % Note that here we adopted the proposed values in Turner et al., 2013,
        % which might partly depend on Gelman et al., 1996,
        % with some adaptation
        gamma=2.38/sqrt(2*Nparam); % scales the difference vector
        eps=0.01; % random jitter (the range of the uniform distribution)
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
        eps=0.01; % scales the constant that avoids the covariance matrix from being singular
    else
        Sd=config.MCMCparam(1);
        t0=config.MCMCparam(2);
        eps=config.MCMCparam(3);
    end
else
    warning('The algorithm is invalid, use DE-MCMC instead.')
    config.Algorithm='DE';
end
% pre-allocation
TrueSamples=zeros(Nchain, Nparam, Nstate);
Samples=zeros(Nchain,Nparam,Nstate);
logPosterior=zeros(Nchain,Nstate);
logPrior=zeros(Nchain,Nstate);
logLikelihood=zeros(Nchain,Nstate);
ConvergenceStat=zeros(1,MAXbatchburnin);
% shuffle random seed
rng(now);

%% MCMC
if strcmp(config.Algorithm,'DE') % Differential Evolution
    if strcmp(config.Verbosity,'iter')
        fprintf('\n\nNow start MCMC sampling based on DE \n\n')
    end
    % set initial start value
    if Transform==1
        start=MCMCConvert_BMW(model.Constraints.start,model.Constraints.ub,model.Constraints.lb,config.Transform); % do transform
    else
        start=model.Constraints.start; % use true parameter value (constrain the transition kernel instead)
    end
    count=0; % number of finished burn-in batches
    while 1
        if count>=MAXbatchburnin
            if strcmp(config.Verbosity,'notify')
                fprintf('\nReached the max number of burn-in samples without convergence...\n\r')
            end
            if strcmp(config.Verbosity,'iter')
                fprintf('\nReached the max number of burn-in samples without convergence...\n\r')
                fprintf('Now start generate valid samples...\n')
            end
            Summary.ConvergenceNbatch=0;
            break;
        end
        % refresh storage
        BurninSamples=zeros(Nchain, Nparam, Nbatchburnin);
        BurninPosterior=zeros(Nchain, Nbatchburnin);
        % generate burnin samples till convergence
        for state=1:Nbatchburnin % loop states
            if state==1
                PrevState=start; % start
            else
                PrevState=BurninSamples(:,:,state-1); % previous state
            end
            for chain=1:Nchain % loop chains
                % generate proposal state
                ChainRange=1:Nchain; % decision space of chains
                ChainRange=ChainRange(ChainRange~=chain); % exclude the current chain
                while 1 % redo sampling if the current sample is out of bounds
                    ChainLot=randsample(ChainRange,2); % ramdomly pick two chains without repetition
                    ChainDiff=gamma*(PrevState(ChainLot(1),:)-PrevState(ChainLot(2),:)); % get difference vector
                    DiffNum=randsample(1:Nparam,1);
                    ChainDiff(DiffNum)=randsample([1,1,1,0],1)*ChainDiff(DiffNum);
                    UniNoise=2*eps.*rand(1,length(eps))-abs(eps); % noise
                    PropState=PrevState(chain,:)+ChainDiff+UniNoise; % proposal state
                    if Transform==0
                        if ~any(PropState>model.Constraints.ub) && ~any(PropState<model.Constraints.lb) % test bounds
                            break;
                        end
                    else
                        break;
                    end
                end
                % test proposal state
                minMHr=rand; % sample Metropolis-Hastings ratio threshold
                switch Transform
                    case 1
                        TruePropState=MCMCConvert_BMW(PropState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                        TruePrevState=MCMCConvert_BMW(PrevState(chain,:),model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                    case 0
                        TruePropState=PropState;
                        TruePrevState=PrevState(chain,:);
                end
                eval(['Prop0=',model.Model,'(TruePropState, data, model);'])
                eval(['Prev0=',model.Model,'(TruePrevState, data, model);'])
                eval(['PropPosterior=-Prop0.',RealOutput,';'])
                eval(['PrevPosterior=-Prev0.',RealOutput,';'])
                MHr=exp(PropPosterior-PrevPosterior); % current MH ratio
                if MHr>minMHr
                    CurrentState=PropState; % accept proposal
                    BurninPosterior(chain,state)=PropPosterior; % record Posterior
                else
                    CurrentState=PrevState(chain,:); % reject, use the previous state instead
                    BurninPosterior(chain,state)=PrevPosterior;
                end
                BurninSamples(chain,:,state)=CurrentState; % record sample value
            end
        end
        count=count+1;
        % diagnose convergence
        [Convergence, Stat]=TestConvergence(BurninSamples,config.Convergence.Tol,config.Convergence.Diagnostic); 
        ConvergenceStat(count)=Stat;
        if strcmp(config.Verbosity,'iter')
            fprintf('%d samples generated, R=%d\n',count*Nburnin, Stat)
        end
        switch Convergence
            case 1
                if strcmp(config.Verbosity,'iter')
                    fprintf('Converged! Now start generate valid samples...\n')
                end
                Summary.ConvergenceNbatch=count; % Record the index of the current batch
                break; % converged! break loop & start collecting valid samples
            case 0 % redo MCMC
                start=BurninSamples(:, :, end);
        end
    end
    Ntrial=length(Prop0.LPPD); % record # trial
    logPointwiseLH=zeros(Nchain,Ntrial,Nstate);
    % generate valid samples after convergence
    restart=BurninSamples(:, :, end);
    for state=1:Nstate % loop states
        for chain=1:Nchain % loop chains
            % generate proposal state
            if state==1
                PrevState=restart; % start
            else
                PrevState=Samples(:,:,state-1); % previous state
            end
            ChainRange=1:Nchain; % decision space of chains
            ChainRange=ChainRange(ChainRange~=chain); % exclude the current chain
            while 1 % redo sampling if the current sample is out of bounds
                ChainLot=randsample(ChainRange,2); % ramdomly pick two chains without repetition
                ChainDiff=gamma*(PrevState(ChainLot(1),:)-PrevState(ChainLot(2),:)); % get difference vector
                DiffNum=randsample(1:Nparam,1);
                ChainDiff(DiffNum)=randsample([0,1,1,1],1)*ChainDiff(DiffNum);
                UniNoise=2*eps.*rand(1,length(eps))-abs(eps); % noise
                PropState=PrevState(chain,:)+ChainDiff+UniNoise; % proposal state
                if Transform==0
                    if ~any(PropState>model.Constraints.ub) && ~any(PropState<model.Constraints.lb) % test bounds
                        break;
                    end
                else
                    break;
                end
            end
            % test proposal state
            switch Transform
                case 1
                    TruePropState=MCMCConvert_BMW(PropState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                    TruePrevState=MCMCConvert_BMW(PrevState(chain,:),model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                case 0
                    TruePropState=PropState;
                    TruePrevState=PrevState(chain,:);
            end
            minMHr=rand; % sample Metropolis-Hastings ratio threshold
            eval(['Prop0=',model.Model,'(TruePropState, data, model);'])
            eval(['Prev0=',model.Model,'(TruePrevState, data, model);'])
            eval(['PropPosterior=-Prop0.',RealOutput,';'])
            eval(['PrevPosterior=-Prev0.',RealOutput,';'])
            MHr=exp(PropPosterior-PrevPosterior); % current MH ratio
            if MHr>minMHr
                CurrentState=PropState; % accept proposal
                if Transform==1, TrueCurrentState=TruePropState; end
                logPosterior(chain,state)=PropPosterior; % record log posterior
                logLikelihood(chain,state)=-Prop0.LLH; % record log likelihood
                logPrior(chain,state)=-Prop0.Prior; % record prior
                logPointwiseLH(chain,:,state)=-Prop0.LPPD; % record pointwise likelihood
            else
                CurrentState=PrevState(chain,:); % reject, use the previous state instead
                if Transform==1, TrueCurrentState=TruePrevState; end
                logPosterior(chain,state)=PrevPosterior; % record log posterior
                logLikelihood(chain,state)=-Prev0.LLH; % record log likelihood
                logPrior(chain,state)=-log(Prev0.Prior); % record prior
                logPointwiseLH(chain,:,state)=-Prev0.LPPD; % record pointwise likelihood
            end
            Samples(chain,:,state)=CurrentState; % record sample value
            if Transform==1, TrueSamples(chain,:,state)=TrueCurrentState; end
        end
        if strcmp(config.Verbosity,'iter') && mod(state*Nchain,500)<Nchain && state*Nchain/500>=1
            fprintf('%d/%d samples collected after convergence...\n', state*Nchain, config.Nsample)
        end
    end
    if strcmp(config.Verbosity,'iter') || strcmp(config.Verbosity,'notify') || strcmp(config.Verbosity,'final')
        fprintf('%d samples collected in total...\n', state*chain)
    end
    
elseif strcmp(config.Algorithm,'MH') % Metropolis-Hastings
    if strcmp(config.Verbosity,'iter')
        fprintf('\n\nNow start MCMC sampling based on MH \n\n')
    end
    % set initial values
    if Transform==1
        start=MCMCConvert_BMW(model.Constraints.start,model.Constraints.ub,model.Constraints.lb,config.Transform); % do transform
    else
        start=model.Constraints.start; % use true parameter value (constrain the transition kernel instead)
    end
    count=0; % number of finished burn-in batches
    while 1
        if count>=MAXbatchburnin
            if strcmp(config.Verbosity,'notify')
                fprintf('Reached the max number of burn-in samples without convergence...\n')
            end
            if strcmp(config.Verbosity,'iter')
                fprintf('Reached the max number of burn-in samples without convergence...\n')
                fprintf('Now start generate valid samples...\n')
            end
            Summary.ConvergenceNbatch=0;
            break;
        end
        BurninSamples=zeros(Nchain, Nparam, Nbatchburnin);
        BurninPosterior=zeros(Nchain, Nbatchburnin);
        % generate burn-in samples till convergence
        state_count=0;
        for state=1:Nbatchburnin
            state_count=state_count+1;
            for chain=1:Nchain
                % update proposal distribution
                if state==1
                    History=start(chain,:); % start
                else
                    % append the past samples to the start value
                    % use inverse to avoid the problem of Nparam and (state-1) being equal
                    History=vertcat(start(chain,:),reshape(BurninSamples(chain,:,1:(state-1)),[Nparam, state-1])');
                end
                if state_count==10 % every ten steps
                    cov=UpdateCov(cov,History,Sd,t0,eps); % adapt the tuning parameters based on the history of sampling
                    state_count=0;
                end
                % generate proposal state
                PrevState=History(end,:); % get previous state
                PropState=mvnrnd(PrevState,cov); % sample from multivariate normal distribution
                % test proposal state
                switch Transform
                    case 1
                        TruePropState=MCMCConvert_BMW(PropState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                        TruePrevState=MCMCConvert_BMW(PrevState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                    case 0
                        TruePropState=PropState;
                        TruePrevState=PrevState;
                end
                minMHr=rand; % sample Metropolis-Hastings ratio threshold
                eval(['Prop0=',model.Model,'(TruePropState, data, model);'])
                eval(['Prev0=',model.Model,'(TruePrevState, data, model);'])
                eval(['PropPosterior=-Prop0.',RealOutput,';'])
                eval(['PrevPosterior=-Prev0.',RealOutput,';'])
                MHr=exp(PropPosterior-PrevPosterior)*(mvnpdf(PrevState,PropState,cov)/mvnpdf(PropState,PrevState,cov)); % current MH ratio
                if MHr>minMHr
                    CurrentState=PropState; % accept proposal
                    BurninPosterior(chain,state)=PropPosterior; % record logPosterior
                else
                    CurrentState=PrevState; % reject, use the previous state instead
                    BurninPosterior(chain,state)=PrevPosterior;
                end
                BurninSamples(chain,:,state)=CurrentState; % record sample value
            end
        end
        count=count+1;
        % diagnose convergence
        [Convergence, Stat]=TestConvergence(BurninSamples,config.Convergence.Tol,config.Convergence.Diagnostic); % use conventional GR statistics here
        if isnan(Stat), Stat=0; end
        ConvergenceStat(count)=Stat;
        if strcmp(config.Verbosity,'iter')
            fprintf('%d burn-in samples generated, R=%d\n',count*Nbatchburnin, Stat)
        end
        switch Convergence
            case 1
                if strcmp(config.Verbosity,'iter')
                    disp('Converged! Now start generate valid samples...\n')
                end
                break; % converged! break loop & start collecting valid samples
            case 0 % refresh storage & redo MCMC
                start=BurninSamples(:, :, end);
        end
    end
    Ntrial=length(Prop0.LPPD); % record # trial
    logPointwiseLH=zeros(Nchain,Ntrial,Nstate);
    % generate valid samples after convergence
    restart=BurninSamples(:, :, end);
    state_count=0;
    for state=1:Nstate
        state_count=state_count+1;
        for chain=1:Nchain
            % update proposal distribution
            if state==1
                History=restart(chain,:); % start
            else
                % append the past samples to the start value
                % use inverse to avoid the problem of Nparam and (state-1) being equal
                History=vertcat(restart(chain,:),reshape(Samples(chain,:,1:(state-1)),[Nparam, state-1])');
            end
            if state_count==10
                cov=UpdateCov(cov,History,Sd,t0,eps); % adapt the tuning parameters based on the history of sampling
                state_count=0;
            end
            % generate proposal state
            PrevState=History(end,:); % get previous state
            PropState=mvnrnd(PrevState,cov); % sample from multivariate normal distribution
            switch Transform
                case 1
                    TruePropState=MCMCConvert_BMW(PropState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                    TruePrevState=MCMCConvert_BMW(PrevState,model.Constraints.ub,model.Constraints.lb,['Inverse', config.Transform]);
                case 0
                    TruePropState=PropState;
                    TruePrevState=PrevState(chain,:);
            end
            % test proposal state
            minMHr=rand; % sample Metropolis-Hastings ratio threshold
            eval(['Prop0=',model.Model,'(TruePropState, data, model);'])
            eval(['Prev0=',model.Model,'(TruePrevState, data, model);'])
            eval(['PropPosterior=-Prop0.',RealOutput,';'])
            eval(['PrevPosterior=-Prev0.',RealOutput,';'])
            MHr=exp(PropPosterior-PrevPosterior)*(mvnpdf(PrevState,PropState,cov)/mvnpdf(PropState,PrevState,cov)); % current MH ratio
            if MHr>minMHr
                CurrentState=PropState; % accept proposal
                if Transform==1, TrueCurrentState=TruePropState; end
                logPosterior(chain,state)=PropPosterior; % record logPosterior
                logLikelihood(chain,state)=-Prop0.LLH; % record likelihood
                logPrior(chain,state)=-Prop0.Prior; % record prior;
                logPointwiseLH(chain,:,state)=-Prop0.LPPD; % record pointwise likelihood
            else
                CurrentState=PrevState; % reject, use the previous state instead
                if Transform==1, TrueCurrentState=TruePrevState; end
                logPosterior(chain,state)=PrevPosterior; % record logPosterior
                logLikelihood(chain,state)=-Prev0.LLH; % record likelihood
                logPrior(chain,state)=-log(Prev0.Prior); % record prior;
                logPointwiseLH(chain,:,state)=-Prev0.LPPD; % record pointwise likelihood
            end
            Samples(chain,:,state)=CurrentState; % record sample value
            if Transform==1, TrueSamples(chain,:,state)=TrueCurrentState; end
        end
        if strcmp(config.Verbosity,'iter') && mod(state*Nchain,500)<Nchain && state*Nchain/500>=1
            fprintf('%d/%d samples collected...\n', state*Nchain, config.Nsample)
        end
    end
    if strcmp(config.Verbosity,'iter') || strcmp(config.Verbosity,'notify') || strcmp(config.Verbosity,'final')
        fprintf('%d samples collected in total...\n', state*chain)
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
[Summary.MAXposterior, IndMAXposterior]=max(RawSampling.logPosterior); % max of logPosterior density
Summary.FitParam=RawSampling.Samples(IndMAXposterior,:); % best parameter(s)

%% Epilogue
if ~strcmp(config.Verbosity,'off')
    fprintf('\nAll done!\n\n')
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
elseif strcmp(method,'GRL'), [BoolConverge, Stat]=TestConvergence_GRL(History, config);...
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
if R>1+delta
    BoolConverge=1; % converged
else
    BoolConverge=0; % not coverged yet
end
end

% TestConvergence_GRL based on the adapted GR statistic in Vats & Knudson, 2019
function [BoolConverge, R]=TestConvergence_GRL(History, config)
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
if R>1+delta
    BoolConverge=1; % converged
else
    BoolConverge=0; % not coverged yet
end
end
