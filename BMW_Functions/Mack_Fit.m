%% Mack Fit
% Fit designated model by MLE/MAP
% Assess designated model by calculating LP
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 8/31/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Param, Quality]=Mack_Fit(Data, Config, Model, Constraints, FitOptions)

if nargin==4 || isempty(FitOptions)
    FitOptions.Algorithm='DE-MCMC'; % set MCMC as default
elseif nargin==3 || isempty(Constraints)
    error('Constraints are needed...')
elseif nargin==2 || isempty(Model)
    error('Cannot find likelihood function...')
end

if strcmp(FitOptions.Algorithm,'fmincon: sqp')
        %% fmincon: sqp
        % Matlab Optimization Toolbox
        % https://ww2.mathworks.cn/help/optim/ug/fmincon.html;jsessionid=e52720ed65cadadf60d137ef2b4c#d117e83832
        if ~exist('fmincon','file')
            error('Error: Optimization toolbox is needed.')
        end
        if ~isfield(FitOptions,'fminconOptions')
            % Algorithm: 'active-set'/'interior-point'/'sqp'/'sqp-legacy'/'trust-region-reflective'
            FitOptions.fminconOptions.Display='iter';
            FitOptions.fminconOptions.MaxIter=5000;
            FitOptions.fminconOptions.StepTolerance=1e-6;
        end
        FitOptions.fminconOptions.Algorithm='sqp';
        eval(['[Param_Mack, LP_Mack, Exitflag, Output]=fmincon(@(Param)',Model,'(Param, Data, Config), Constraints.start,[],[],[],[],Constraints.lb, Constraints.ub, [], FitOptions.fminconOptions);'])

elseif strcmp(FitOptions.Algorithm,'fmincon: interior-point')
    %% fmincon: interior-point
        % Matlab Optimization Toolbox
        % https://ww2.mathworks.cn/help/optim/ug/fmincon.html;jsessionid=e52720ed65cadadf60d137ef2b4c#d117e83832
        if ~exist('fmincon','file')
            error('Error: Optimization toolbox is needed.')
        end
        if ~isfield(FitOptions,'fminconOptions')
            % Algorithm: 'active-set'/'interior-point'/'sqp'/'sqp-legacy'/'trust-region-reflective'
            FitOptions.fminconOptions.Display='iter';
            FitOptions.fminconOptions.MaxIter=3000;
            FitOptions.fminconOptions.StepTolerance=1e-6;
        end
        FitOptions.fminconOptions.Algorithm='interior-point';
        eval(['[Param_Mack, LP_Mack, Exitflag, Output]=fmincon(@(Param)',Model,'(Param, Data, Config), Constraints.start,[],[],[],[],Constraints.lb, Constraints.ub, [], FitOptions.fminconOptions);'])

elseif strcmp(FitOptions.Algorithm,'fmincon: active-set')
    %% fmincon: active-set
        % Matlab Optimization Toolbox
        % https://ww2.mathworks.cn/help/optim/ug/fmincon.html;jsessionid=e52720ed65cadadf60d137ef2b4c#d117e83832
        if ~exist('fmincon','file')
            error('Error: Optimization toolbox is needed.')
        end
        if ~isfield(FitOptions,'fminconOptions')
            % Algorithm: 'active-set'/'interior-point'/'sqp'/'sqp-legacy'/'trust-region-reflective'
            FitOptions.fminconOptions.Display='iter';
            FitOptions.fminconOptions.MaxIter=3000;
            FitOptions.fminconOptions.StepTolerance=1e-6;
        end
        FitOptions.fminconOptions.Algorithm='active-set';
        eval(['[Param_Mack, LP_Mack, Exitflag, Output]=fmincon(@(Param)',Model,'(Param, Data, Config), Constraints.start,[],[],[],[],Constraints.lb, Constraints.ub, [], FitOptions.fminconOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'BADS')
        %% bads
        % Acerbi & Ma, 2017, Advances in Neural Information Processing Systems
        % http://github.com/lacerbi/bads
        if ~exist('bads','file')
            error('Error: BADS toolbox is needed.')
        end
        if ~isfield(FitOptions,'badsOptions')
            FitOptions.badsOptions=bads('defaults');
        end
        eval(['[Param_Mack, LP_Mack, Exitflag, Output] = bads(@(Param)' Model, '(Param, Data, Config), Constraints.start, Constraints.lb, Constraints.ub, Constraints.lb, Constraints.ub, [], FitOptions.badsOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'MADS')
        %% mads
        % Audet & Dennis, 2006, SIAM Journal on Optimization
        % Matlab Global Optimization Toolbox
        % https://ww2.mathworks.cn/help/gads/index.html?s_tid=CRUX_lftnav
        if ~exist('patternsearch','file')
            error('Error: Global Optimization toolbox is needed.')
        end
        eval(['[Param_Mack, LP_Mack, Exitflag, Output] = patternsearch(@(Param)' Model, '(Param, Data, Config), Constraints.start, [], [], [], [], Constraints.lb, Constraints.ub, [], FitOptions.madsOptions);'])
                
elseif strcmp(FitOptions.Algorithm,'GA')
        %% Genetic Algorithm
        % Matlab Global Optimization Toolbox
        % https://ww2.mathworks.cn/help/gads/index.html?s_tid=CRUX_lftnav
        if ~exist('ga','file')
            error('Error: Global Optimization toolbox is needed.')
        end
        eval(['[Param_Mack, LP_Mack, Exitflag, Output] = ga(@(Param)' Model, '(Param, Data, Config), length(Constraints.start), [], [], [], [], Constraints.lb, Constraints.ub, [], FitOptions.gaOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'SA')
        %% Simulated Annealing
        % Matlab Global Optimization Toolbox
        % https://ww2.mathworks.cn/help/gads/index.html?s_tid=CRUX_lftnav
        if ~exist('simulannealbnd','file')
            error('Error: Global Optimization toolbox is needed.')
        end
        eval(['[Param_Mack, LP_Mack, Exitflag, Output] = simulannealbnd(@(Param)' Model, '(Param, Data, Config), Constraints.start, Constraints.lb, Constraints.ub, FitOptions.saOptions);'])
        
elseif strcmp(FitOptions.Algorithm,'DE-MCMC')
        %% DE-MCMC
        % Default algorithm
        % Differential Evolution Monte Carlo Markov Chain
        % Built-in function in Bayesian Modeling of Working Memory (BMW) Toolbox
        if ~exist('Mack_MCMC','file')
            error('Error: Mack_MCMC function not detected.')
        end
        Model_MCMC=Config;
        Model_MCMC.Model=Model;
        Model_MCMC.Constraints=Constraints;
        % check start values
        if size(Model_MCMC.Constraints.start,1)==1
            Model_MCMC.Constraints.start=repmat(Model_MCMC.Constraints.start,[4,1]);
        end
        Config_MCMC.Algorithm='DE';
        if isfield(FitOptions,'MCMCoptions')
            Config_MCMC=FitOptions.MCMCoptions;
        end
        [MCMCResult,OptResult]=Mack_MCMC(Model_MCMC, Data, Config_MCMC);
        Param_Mack=OptResult.FitParam;
        LP_Mack=log(OptResult.MAXposterior);
        Quality.MCMCResult=MCMCResult;
elseif strcmp(FitOptions.Algorithm,'MH-MCMC')
        %% MH-MCMC
        % (Adaptive) Metropolis-Hastings Monte Carlo Markov Chain
        % Built-in function in Bayesian Modeling of Working Memory (BMW) Toolbox
        if ~exist('Mack_MCMC','file')
            error('Error: Mack_MCMC function not detected.')
        end
        Model_MCMC=Config;
        Model_MCMC.Model=Model;
        Model_MCMC.Constraints=Constraints;
        % check start values
        if size(Model_MCMC.Constraints.start,1)==1
            Model_MCMC.Constraints.start=repmat(Model_MCMC.Constraints.start,[2*size(Model_MCMC.Constraints.start,2)+1,1]);
        end
        Config_MCMC.Algorithm='MH';
        if isfield(FitOptions,'MCMCoptions')
            Config_MCMC=FitOptions.MCMCoptions;
        end
        [MCMCResult,OptResult]=Mack_MCMC(Model_MCMC, Data, Config_MCMC);
        Param_Mack=OptResult.FitParam;
        LP_Mack=log(OptResult.MAXposterior);
        Quality.MCMCResult=MCMCResult;
end

if ~exist('LP_Mack','var')
    error('Sorry, the algorithm is invalid...')
end

Quality.LP=-LP_Mack;
Param=Param_Mack;

end
