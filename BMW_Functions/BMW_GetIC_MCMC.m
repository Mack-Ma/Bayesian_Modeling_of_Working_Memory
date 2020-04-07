%% Get ICs based on posterior samples
%
% Extract the information criterion of interest from the MCMC outputs
%
% -----------------------
% IC=BMW_GetIC_MCMC(Samples, Model, Data, Method)
%
% ## Input ##
% - Samples
%   mat, contains posterior samples
%   help BMW_MCMC for details
%
% - Model
%   struct, configures the model
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
%
% - Data
%   mat, the data variable that qualifies the simulation.
%   Type BMW('manual') for detailed requirements.
%
% - Method
%   string, designates the method to compare models
%       'DIC', Deviance Information Criterion
%       'DIC*', modified Deviance Information Criterion, which gives a larger
%       penalty on model complexity than DIC
%       'WAIC1', Watanabe-Akaike Information Criterion, which uses the 
%       difference between pointwise mean deviance and pointwise deviance 
%       of mean as penalty
%       'WAIC2', Watanabe-Akaike Information Criterion, which uses the
%       sum of the pointwise variance of likelihood as penalty
%       'LME_HarmonicMean', log model evidence, calculated by the
%       generalized harmonic mean estimator
%       'LME_BridgeSampling', log model evidence, calculated by the
%       bridge sampling estimator
%
% ## Output ##
% - IC
% float, return the information criterion of interest
%
% ## Reference ##
% - 
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 3/22/2020
%
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/BMW-Ma/Bayesian_Modeling_of_Working_Memory
%

function IC=BMW_GetIC_MCMC(RawSampling, Model, Data, Method)
RawSamples=RawSampling.RawSamples;
Samples=RawSampling.Samples;
Posterior=RawSampling.logPosterior;
Likelihood=RawSampling.logLikelihood;
LPPD=RawSampling.logPointwiseLikelihood;
Nsample=size(Samples,1);
Nparam=size(Samples,2);
% set default
if nargin==3 || ~isfield(Method,'IC')
    Method.IC='LME_HarmonicMean';
end
if nargin==3 || ~isfield(Method, 'Verbosity')
    Method.Verbosity=1;
end
% get IC
switch Method.IC
    case 'DIC' % get deviance information criterion
        if Method.Verbosity==1;
            fprintf('\nNow start estimating DIC... \n\n')
        end
        % get likelihood
        LLH=-Likelihood;
        % get posterior mean deviance (expected log likelihood)
        MDev=mean(2*LLH);
        % get deviance of mean
        MeanSample=mean(Samples,1);
        eval(['DevM=2*',Model.Model,'(MeanSample,Data,Model);'])
        % get 'effective number of parameters'
        pD=MDev-DevM;
        % DIC
        IC=MDev+pD;
        if Method.Verbosity==1;
            fprintf('Done!\n')
        end
    case 'DIC*' % the only difference between DIC & DIC* is that DIC* gives larger penalty to model complexity
        if Method.Verbosity==1;
            fprintf('\nNow start estimating DIC*... \n\n')
        end
        % get likelihood
        LLH=-Likelihood;
        % get posterior mean deviance (expected log likelihood)
        MDev=mean(2*LLH);
        % get deviance of mean
        MeanSample=mean(Samples,1);
        eval(['DevM=2*',Model.Model,'(MeanSample,Data,Model);'])
        % get 'effective number of parameters'
        pD=MDev-DevM;
        % DIC*
        IC=MDev+2*pD;
        if Method.Verbosity==1;
            fprintf('Done!\n')
        end
    case 'WAIC1' % get watanabe-akaike information criterion
        if Method.Verbosity==1;
            fprintf('\nNow start estimating WAIC1... \n\n')
        end
        % get log pointwise predictive density
        lppd=LPPD;
        elppd=-2*sum(log(mean(exp(lppd),1)));
        % get penalty
        p_waic1=2*sum(log(mean(exp(lppd),1)))-2*sum(mean(lppd,1));
        % WAIC
        IC=elppd-p_waic1;
        if Method.Verbosity==1;
            fprintf('Done!\n')
        end
    case 'WAIC2' % get watanabe-akaike information criterion
        if Method.Verbosity==1;
            fprintf('\nNow start estimating WAIC2... \n\n')
        end
        % get log pointwise predictive density
        lppd=LPPD;
        elppd=-2*sum(log(mean(exp(lppd),1)));
        % get penalty
        p_waic2=0.5*sum(var(lppd,1));
        % WAIC
        IC=elppd-p_waic2;
        if Method.Verbosity==1;
            fprintf('Done!\n')
        end
    case 'LME_HarmonicMean' % get log model evidence (marginal likelihood) through the generalized harmonic mean estimator
        if Method.Verbosity==1;
            fprintf('\nNow estimate marginal likelihood based on the generalized harmonic mean estimator... \n\n')
        end
        % I hate the dogmaticity that the importance density brought in.
        % Nonetheless, it's definitely faster than bridge sampling.
        % use the ML estimator to find the fittest proposal distribution
        % get the covariance matrix of the raw posterior samples
        SampleCov=cov(RawSamples);
        IDMean=mean(RawSamples);
        IDCov=SampleCov*(Nsample-1)/Nsample; % Unbiased estimator
        % get log importance density
        LP=Posterior;
        LID=zeros(Nsample,1);
        for i=1:Nsample
            LID(i)=log(mvnpdf(RawSamples(i,:),IDMean,IDCov)); % we use multivariate normal distribution
        end
        % scaling
        LID=LID*(max(LP)-500/max(LID));
        % ID/posterior
        DR=exp(LID-LP);
        % test bounds (due to the computational limit)
        DR(DR==Inf)=realmax('double')/Nsample;
        % harmonic mean estimator
        IC=-log((mean(DR)));
        if Method.Verbosity==1
            fprintf('Done!\n')
        end
    case 'LME_BridgeSampling' % get log model evidence (marginal likelihood) through bridge sampling
        if Method.Verbosity==1;
            fprintf('\nNow estimate marginal likelihood based on bridge sampling... \n\n')
        end
        % recommended method, takes time tho
        % hate the dogmaticity as well....
        % use ML estimator to fit the proposal distribution (mvn)
        if Method.Verbosity==1;
            fprintf('Constructing the proposal function... \n')
        end
        % get the covariance matrix of the raw posterior samples
        SampleCov=cov(RawSamples);
        IDCov=SampleCov*(Nsample-1)/Nsample; % Unbiased estimator
        % Sampling based on the proposal distribution
        Nprop=Nsample/2;
        [PropSamples,PropID]=randmvn(Nprop,zeros(1,Nparam),IDCov);
        if Method.Verbosity==1;
            fprintf('Done!\n')
            fprintf('Now start updating the estimated marginal likelihood value... \n')
        end
        % get log posterior of the proposal samples
        % do transform
        TruePropSamples=MCMCConvert_BMW(PropSamples,Model.Constraints.ub,Model.Constraints.lb,'InverseProbit');
        LPprop=zeros(Nprop,1);
        Model.Output='LP';
        for i=1:Nprop
            eval(['LPprop(i)=-',Model.Model,'(TruePropSamples(i,:), Data, Model);'])
        end
        % get proposal density of the target samples
        TarID=mvnpdf(RawSamples,[],IDCov);
        % get L1 & L2
        LL1=Posterior-log(TarID);
        LL2=LPprop-log(PropID);
        LLstar=median(LL1); % avoid reaching the computational limit, see https://osf.io/8x7m9/
        % use the optimal bridging function in Gronau et al, 2017
        MaxT=10000;
        ToleranceME=1e-10;
        ME0=0; % initial guess
        s1=Nsample/(Nsample+Nprop);
        s2=Nprop/(Nsample+Nprop);
        for t=1:MaxT
            ME1=sum(exp(LL2-LLstar)./(s1*exp(LL2-LLstar)+s2*ME0));
            ME2=sum(1./(s1*exp(LL1-LLstar)+s2*ME0));
            if ME2==Inf
                ME=realmin('double');
            else
                ME=(s1/s2)*ME1/ME2;
            end
            % check tolerance
            if abs(1-ME0/ME)<=ToleranceME
                MEfinal=ME;
                fprintf('Converged!\n')
                break;
            end
            ME0=ME;
        end
        if t==MaxT
            fprintf('Reach max number of interation, marginal likelihood didn''t converge...\n')
            MEfinal=ME;
        end
        IC=log(MEfinal)+LLstar;
end

end

% draw samples from mvn distribution based on DE-MCMC
function [AllSamples, AllDensity]=randmvn(N,mean,covar)
Nchain=5;
Nsample=ceil(N/Nchain);
Samples=zeros(Nchain,length(mean),Nsample);
Density=zeros(Nchain,Nsample);
% construct initial values
start=zeros(Nchain,length(mean));
for chain=1:Nchain
    start(chain,:)=mean+rand(1,length(mean));
end
% DE parameters
gamma=2.38/sqrt(length(mean));
eps=0.001;
for state=1:Nsample
    if state==1
        PrevState=start;
    else
        PrevState=Samples(:, :, state-1);
    end
    for chain=1:Nchain
        % generate proposal state
        ChainRange=1:Nchain;
        ChainRange(chain)=[];
        ChainInd=randsample(ChainRange,2);
        ChainDiff=gamma*(PrevState(ChainInd(1),:)-PrevState(ChainInd(2),:)); % randomly choose two chains
        UniformNoise=2*eps*rand(1,length(mean))-eps; % noise
        PropState=PrevState(chain,:)+ChainDiff+UniformNoise;
        % test proposal state
        minMHr=rand; % generate MH ratio threshold
        PropDensity=mvnpdf(PropState,mean,covar);
        PrevDensity=mvnpdf(PrevState(chain,:),mean,covar);
        curMHr=PropDensity/PrevDensity; % current MH ratio
        if curMHr>minMHr
            Samples(chain,:,state)=PropState;
            Density(chain,state)=PropDensity;
        else
            Samples(chain,:,state)=PrevState(chain,:);
            Density(chain,state)=PrevDensity;
        end
    end
end

% merge chains
AllSamples=zeros(Nsample*Nchain,length(mean));
AllDensity=zeros(Nsample*Nchain,1);
for chain=1:Nchain
    AllDensity((1+(chain-1)*Nsample):chain*Nsample)=Density(chain,:)';
    AllSamples((1+(chain-1)*Nsample):chain*Nsample,:)=permute(Samples(chain,:,:),[3, 2, 1]);
end
AllDensity=AllDensity(end-N+1:end,:);
AllSamples=AllSamples(end-N+1:end,:);
end