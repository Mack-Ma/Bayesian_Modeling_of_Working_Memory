%% BMW Fit
% Fit designated model by MLE/MAP
% Assess designated model by calculating LP
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 8/31/2019
%
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox: https://github.com/BMW-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Param, Quality]=BMW_Fit(Data, Config, Model, Constraints, FitOptions)

rng(now); % refresh the random seed
if nargin==4 || isempty(FitOptions) || ~isfield(FitOptions,'Algorithm')
    FitOptions.Algorithm='DE-MCMC'; % set MCMC as default
elseif nargin==3 || isempty(Constraints)
    error('Constraints are needed...')
elseif nargin==2 || isempty(Model)
    error('Cannot detect the likelihood function...')
end

if strcmp(FitOptions.Algorithm,'fmincon: sqp')
        %% fmincon: sqp
        % Matlab Optimization Toolbox
        % https://ww2.mathworks.cn/help/optim/ug/fmincon.html;jsessionid=e52720ed65cadadf60d137ef2b4c#d117e83832
        if ~exist('fmincon','file')
            error('Error: Optimization toolbox is needed.')
        end
        FitOptions.fminconOptions.Verbosity=FitOptions.Display;
        FitOptions.fminconOptions.Algorithm='sqp';
        eval(['[Param_BMW, LP_BMW, Exitflag, Output]=fmincon(@(Param)',Model,'(Param, Data, Config), Constraints.start,[],[],[],[],Constraints.lb, Constraints.ub, [], FitOptions.fminconOptions);'])

elseif strcmp(FitOptions.Algorithm,'fmincon: interior-point')
    %% fmincon: interior-point
        % Matlab Optimization Toolbox
        % https://ww2.mathworks.cn/help/optim/ug/fmincon.html;jsessionid=e52720ed65cadadf60d137ef2b4c#d117e83832
        if ~exist('fmincon','file')
            error('Error: Optimization toolbox is needed.')
        end
        FitOptions.fminconOptions.Verbosity=FitOptions.Display;
        FitOptions.fminconOptions.Algorithm='interior-point';
        eval(['[Param_BMW, LP_BMW, Exitflag, Output]=fmincon(@(Param)',Model,'(Param, Data, Config), Constraints.start,[],[],[],[],Constraints.lb, Constraints.ub, [], FitOptions.fminconOptions);'])

elseif strcmp(FitOptions.Algorithm,'fmincon: active-set')
    %% fmincon: active-set
        % Matlab Optimization Toolbox
        % https://ww2.mathworks.cn/help/optim/ug/fmincon.html;jsessionid=e52720ed65cadadf60d137ef2b4c#d117e83832
        if ~exist('fmincon','file')
            error('Error: Optimization toolbox is needed.')
        end
        FitOptions.fminconOptions.Verbosity=FitOptions.Display;
        FitOptions.fminconOptions.Algorithm='active-set';
        eval(['[Param_BMW, LP_BMW, Exitflag, Output]=fmincon(@(Param)',Model,'(Param, Data, Config), Constraints.start,[],[],[],[],Constraints.lb, Constraints.ub, [], FitOptions.fminconOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'BADS')
        %% bads
        % Acerbi & Ma, 2017, Advances in Neural Information Processing Systems
        % http://github.com/lacerbi/bads
        if ~exist('bads','file')
            error('Error: BADS toolbox is needed.')
        end
        if ~isfield(FitOptions,'BADSOptions')
            FitOptions.BADSOptions=bads('defaults');
        end
        eval(['[Param_BMW, LP_BMW, Exitflag, Output] = bads(@(Param)' Model, '(Param, Data, Config), Constraints.start, Constraints.lb, Constraints.ub, Constraints.lb, Constraints.ub, [], FitOptions.BADSOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'MADS')
        %% mads
        % Audet & Dennis, 2006, SIAM Journal on Optimization
        % Matlab Global Optimization Toolbox
        % https://ww2.mathworks.cn/help/gads/index.html?s_tid=CRUX_lftnav
        if ~exist('patternsearch','file')
            error('Error: Global Optimization toolbox is needed.')
        end
        FitOptions.MADSOptions.Verbosity=FitOptions.Display;
        eval(['[Param_BMW, LP_BMW, Exitflag, Output] = patternsearch(@(Param)' Model, '(Param, Data, Config), Constraints.start, [], [], [], [], Constraints.lb, Constraints.ub, [], FitOptions.MADSOptions);'])
                
elseif strcmp(FitOptions.Algorithm,'GA')
        %% Genetic Algorithm
        % Matlab Global Optimization Toolbox
        % https://ww2.mathworks.cn/help/gads/index.html?s_tid=CRUX_lftnav
        if ~exist('ga','file')
            error('Error: Global Optimization toolbox is needed.')
        end
        FitOptions.GAOptions.Verbosity=FitOptions.Display;
        eval(['[Param_BMW, LP_BMW, Exitflag, Output] = ga(@(Param)' Model, '(Param, Data, Config), length(Constraints.start), [], [], [], [], Constraints.lb, Constraints.ub, [], FitOptions.GAOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'SA')
        %% Simulated Annealing
        % Matlab Global Optimization Toolbox
        % https://ww2.mathworks.cn/help/gads/index.html?s_tid=CRUX_lftnav
        if ~exist('simulannealbnd','file')
            error('Error: Global Optimization toolbox is needed.')
        end
        FitOptions.SAOptions.Verbosity=FitOptions.Display;
        eval(['[Param_BMW, LP_BMW, Exitflag, Output] = simulannealbnd(@(Param)' Model, '(Param, Data, Config), Constraints.start, Constraints.lb, Constraints.ub, FitOptions.SAOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'DE-MCMC')
        %% DE-MCMC
        % Default algorithm
        % Differential Evolution Monte Carlo Markov Chain
        % Built-in function in Bayesian Modeling of Working Memory (BMW) Toolbox
        if ~exist('BMW_MCMC','file')
            error('Error: BMW_MCMC function not detected.')
        end
        Model_MCMC=Config;
        Model_MCMC.Model=Model;
        Model_MCMC.Constraints=Constraints;
        % check start values
        if size(Model_MCMC.Constraints.start,1)<=2
            start_update=zeros(4,size(Model_MCMC.Constraints.start,2));
            start_update(1:size(Model_MCMC.Constraints.start,1),:)=Model_MCMC.Constraints.start;
            for chain=1:4-size(Model_MCMC.Constraints.start,1)
                start_update(size(Model_MCMC.Constraints.start,1)+chain,:)=...
                    Model_MCMC.Constraints.start(randsample(1:size(Model_MCMC.Constraints.start,1),1),:)+...
                    0.02*rand(1,size(Model_MCMC.Constraints.start,2))-0.01;
            end
            Model_MCMC.Constraints.start=start_update;
        end
        if isfield(FitOptions,'MCMCOptions')
            Config_MCMC=FitOptions.MCMCOptions;
        end
        Config_MCMC.Algorithm='DE';
        [MCMCResult,OptResult]=BMW_parMCMC(Model_MCMC, Data, Config_MCMC);
        Param_BMW=OptResult.FitParam;
        LP_BMW=OptResult.BestPosterior;
        Quality.MCMCResult=MCMCResult;
elseif strcmp(FitOptions.Algorithm,'MH-MCMC')
        %% MH-MCMC
        % (Adaptive) Metropolis-Hastings Monte Carlo Markov Chain
        % Built-in function in Bayesian Modeling of Working Memory (BMW) Toolbox
        if ~exist('BMW_MCMC','file')
            error('Error: BMW_parMCMC function not detected.')
        end
        Model_MCMC=Config;
        Model_MCMC.Model=Model;
        Model_MCMC.Constraints=Constraints;
        % check start values
        if size(Model_MCMC.Constraints.start,1)<=2
            start_update=zeros(4,size(Model_MCMC.Constraints.start,2));
            start_update(1:size(Model_MCMC.Constraints.start,1),:)=Model_MCMC.Constraints.start;
            for chain=1:4-size(Model_MCMC.Constraints.start,1)
                start_update(size(Model_MCMC.Constraints.start,1)+chain,:)=...
                    Model_MCMC.Constraints.start(randsample(1:size(Model_MCMC.Constraints.start,1),1),:)+...
                    0.02*rand(1,size(Model_MCMC.Constraints.start,2))-0.01;
            end
            Model_MCMC.Constraints.start=start_update;
        end
        if isfield(FitOptions,'MCMCOptions')
            Config_MCMC=FitOptions.MCMCOptions;
        end
        Config_MCMC.Algorithm='MH';
        [MCMCResult,OptResult]=BMW_parMCMC(Model_MCMC, Data, Config_MCMC);
        Param_BMW=OptResult.FitParam;
        LP_BMW=OptResult.BestPosterior;
        Quality.MCMCResult=MCMCResult;
end

if ~exist('LP_BMW','var')
    error('Sorry, the algorithm is invalid...')
end

Quality.Output=-LP_BMW;
Param=Param_BMW;

end
