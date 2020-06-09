%% ModelFit_BMW
%
% Fit & assess the model(s) designated by the input
% ------------
% MA_BMW=ModelFit_BMW(Input)
%
% For detailed information, plz refer to the manual (type BMW('Manual') in the command line)
% For the complete example, plz refer to "SampleScript1" in the toolbox
% ## Input ##
% - Input
%   cell array, with each element as the output of Configuration_BMW for each
%   model.
%
% ## Output ##
% - MA_BMW
%   struct, with one subfield (MA_BMW.Param) contains the fittest
%   parameters, some subfields contains the configurative options of "Input"
%   and other subfields comprises the results of the model assessment
%   based on each model criterion.
%
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [MA_BMW]=ModelFit_BMW(Input)

fprintf('\n---------------------------------------- \n')
fprintf('\nNow start model fitting... \n')
if ischar(Input)
    load(Input);
elseif isstruct(Input)
    MA_BMW=Input;
else
    error('Sorry, the input file is not valid...')
end
if ~isfield(MA_BMW.FitOptions,'UniformPrior')
    MA_BMW.FitOptions.UniformPrior=1;
end

%% Fitting
Data=MA_BMW.Data;
% Nss=length(unique(Data{1}.SS));
Nsubj=length(Data);
Quality_Fit=cell(Nsubj,1); Output=zeros(Nsubj,1);
LLH=zeros(Nsubj,1); AIC=zeros(Nsubj,1); BIC=zeros(Nsubj,1);
AICc=zeros(Nsubj,1); DIC1=zeros(Nsubj,1); DIC2=zeros(Nsubj,1);
DICs=zeros(Nsubj,1); WAIC1=zeros(Nsubj,1); WAIC2=zeros(Nsubj,1);
LME_GHM=zeros(Nsubj,1); LME_BS=zeros(Nsubj,1);
% Set prior
Config_Fit=MA_BMW.Model;
if strcmp(MA_BMW.FitOptions.Method,'MAP')
    if MA_BMW.FitOptions.UniformPrior==1
        Config_Fit.Output='LLH';
    else
        Config_Fit.Output='LP';
    end
elseif strcmp(MA_BMW.FitOptions.Method,'MLE')
    Config_Fit.Output='LLH';
else
    error('Invalid Fit Method...')
end
% Convert model names
ModelName=MN_Convert_BMW(MA_BMW.Model.Model);
% Fit
fprintf('\nSubject: %s\n','1')
[Param0,Quality]=BMW_Fit(Data{1},Config_Fit,ModelName,MA_BMW.Constraints,MA_BMW.FitOptions);
Ntrial=length(Data{1}.sample);
Nparam=length(Param0);
Param=zeros(Nsubj,Nparam);
Param(1,:)=Param0;
Quality_Fit{1}=Quality;
Output(1)=Quality.Output;
if Nsubj>1
    for subj=2:Nsubj
        fprintf('\nSubject: %s\n',num2str(subj))
        % Fit
        [Param(subj,:),Quality]=BMW_Fit(Data{subj},Config_Fit,ModelName,MA_BMW.Constraints,MA_BMW.FitOptions);
        Quality_Fit{subj}=Quality;
        Output(subj)=Quality.Output;
    end
end
if strcmp(MA_BMW.FitOptions.Method,'LP') && MA_BMW.FitOptions.UniformPrior==1
    Config_Fit.Output='LP';
end

%% LLH, AIC, AICc, BIC
if any(strcmp(MA_BMW.Criteria,'LLH')) || any(strcmp(MA_BMW.Criteria,'AIC')) || any(strcmp(MA_BMW.Criteria,'AICc')) || any(strcmp(MA_BMW.Criteria,'BIC'))
    fprintf('\nNow start estimating LLH/AIC/AICc/BIC...\n')
    % Calculate LLH/AIC/AICc/BIC based on maximum likelihood
    if isfield(Config_Fit,'Output') && (strcmp(Config_Fit.Output,'LLH') || (strcmp(Config_Fit.Output,'LP') && MA_BMW.FitOptions.UniformPrior==1)) % test mode
        for subj=1:Nsubj
            fprintf('\nSubject: %s\n',num2str(subj))
            if any(strcmp(MA_BMW.Criteria,'LLH')), LLH(subj)=-Output(subj); end
            if any(strcmp(MA_BMW.Criteria,'AIC')), AIC(subj)=-2*Output(subj)+2*Nparam; end
            if any(strcmp(MA_BMW.Criteria,'AICc')), AICc(subj)=-2*Output(subj)+2*Nparam+2*Nparam*(Nparam+1)/(Ntrial-Nparam-1); end
            if any(strcmp(MA_BMW.Criteria,'BIC')), BIC(subj)=-2*Output(subj)+Nparam*log(Ntrial); end
        end
    else % find max LLH instead
        fprintf('\nMaximum likelihood is required for LLH/AIC/AICc/BIC...\n')
        MoveOn=input('Do you want to fit model(s) again based on MLE? (y/n)\n','s');
        if strcmp(MoveOn,'y')
            for subj=1:Nsubj
                fprintf('\nSubject: %s\n',num2str(subj))
                % Fit
                [~,Quality_MLE]=BMW_Fit(Data{subj},Config_Fit,ModelName,MA_BMW.Constraints,MA_BMW.FitOptions);
                Output_MLE=Quality_MLE.Output;
                if any(strcmp(MA_BMW.Criteria,'LLH')), LLH(subj)=-Output_MLE; end
                if any(strcmp(MA_BMW.Criteria,'AIC')), AIC(subj)=-2*Output_MLE+2*Nparam; end
                if any(strcmp(MA_BMW.Criteria,'AICc')), AICc(subj)=-2*Output_MLE+2*Nparam+2*Nparam*(Nparam+1)/(Ntrial-Nparam-1); end
                if any(strcmp(MA_BMW.Criteria,'BIC')), BIC(subj)=-2*Output_MLE+Nparam*log(Ntrial); end
            end
        end
    end
    fprintf('\nDone!\n')
end

%% LME (BS/GHM), DIC1, DIC2, DIC*, WAIC1, WAIC2
% Calculate LME based on MCMC
if any(strcmp(MA_BMW.Criteria,'DIC1')) || any(strcmp(MA_BMW.Criteria,'DIC2')) || any(strcmp(MA_BMW.Criteria,'DIC*')) || any(strcmp(MA_BMW.Criteria,'WAIC1'))...
        || any(strcmp(MA_BMW.Criteria,'WAIC2')) || any(strcmp(MA_BMW.Criteria,'LME_BS')) || any(strcmp(MA_BMW.Criteria,'LME_GHM'))
    fprintf('\nNow start estimate LME_BS/LME_GHM/DIC1/DIC2/DIC*/WAIC1/WAIC2...\n')
    if isfield(Quality,'MCMCResult')
        for subj=1:Nsubj
            fprintf('\nSubject: %d',subj)
            Model_MCMC=Config_Fit;
            Model_MCMC.Model=ModelName;
            Model_MCMC.Constraints=MA_BMW.Constraints;
            MCMCcur=Quality_Fit{subj}.MCMCResult;
            switch MA_BMW.FitOptions.Display
                case 'iter'
                    Method.Verbosity=1;
                case 'detail'
                    Method.Verbosity=2;
                case 'off'
                    Method.Verbosity=0;
            end
            if any(strcmp(MA_BMW.Criteria,'LME_BS'))
                Method.IC='LME_BridgeSampling';
                LME_BS(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
            if any(strcmp(MA_BMW.Criteria,'LME_GHM'))
                Method.IC='LME_HarmonicMean';
                LME_GHM(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
            if any(strcmp(MA_BMW.Criteria,'DIC1'))
                Method.IC='DIC1';
                DIC1(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
            if any(strcmp(MA_BMW.Criteria,'DIC2'))
                Method.IC='DIC2';
                DIC2(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
            if any(strcmp(MA_BMW.Criteria,'DIC*'))
                Method.IC='DIC*';
                DICs(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
            if any(strcmp(MA_BMW.Criteria,'WAIC1'))
                Method.IC='WAIC1';
                WAIC1(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
            if any(strcmp(MA_BMW.Criteria,'WAIC2'))
                Method.IC='WAIC2';
                WAIC2(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
            end
        end
    else % do DE-MCMC to collect samples
        fprintf('\nBMW didn''t detect MCMC posterior samples...\n')
        MoveOn=input('\nDo you want redo MCMC? (y/n)\n','s');
        if strcmp(MoveOn,'y')
            for subj=1:Nsubj
                fprintf('Subject: %d\n',subj)
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
                [MCMCcur,~]=BMW_parMCMC(Model_MCMC, Data, Config_MCMC);
                Quality.MCMCResult=MCMCcur;
                switch MA_BMW.FitOptions.Display
                    case 'iter'
                        Method.Verbosity=1;
                    case 'detail'
                        Method.Verbosity=2;
                    case 'off'
                        Method.Verbosity=0;
                end
                if any(strcmp(MA_BMW.Criteria,'LME_GHM'))
                    Method.IC='LME_HarmonicMean';
                    LME_GHM(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if any(strcmp(MA_BMW.Criteria,'LME_BS'))
                    Method.IC='LME_BridgeSampling';
                    LME_BS(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if any(strcmp(MA_BMW.Criteria,'DIC1'))
                    Method.IC='DIC1';
                    DIC1(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if any(strcmp(MA_BMW.Criteria,'DIC2'))
                    Method.IC='DIC2';
                    DIC2(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if any(strcmp(MA_BMW.Criteria,'DIC*'))
                    Method.IC='DIC*';
                    DICs(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if any(strcmp(MA_BMW.Criteria,'WAIC1'))
                    Method.IC='WAIC1';
                    WAIC1(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if any(strcmp(MA_BMW.Criteria,'WAIC2'))
                    Method.IC='WAIC2';
                    WAIC2(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
            end
        end
    end
    fprintf('\nDone!\n')
end

%% Epilogue
MA_BMW.Param=Param;
MCMCAll=cell(1,Nsubj);
if isfield(Quality,'MCMCResult')
    for subj=1:Nsubj
        MCMCAll{subj}=Quality_Fit{subj}.MCMCResult;
    end
end
MA_BMW.MCMC=MCMCAll;
if isfield(MA_BMW.Criteria,'Output') && MA_BMW.Criteria.Output==1
    MA_BMW.Output=Output;
end
C.LLH=LLH;
C.AIC=AIC;
C.AICc=AICc;
C.BIC=BIC;
C.DIC1=DIC1;
C.DIC2=DIC2;
C.DICs=DICs;
C.WAIC1=WAIC1;
C.WAIC2=WAIC2;
C.LME_BS=LME_BS;
C.LME_GHM=LME_GHM;
% delete all empty subfields
MA_MC=fieldnames(C);
for i=1:length(MA_MC)
    eval(['CurField=C.',MA_MC{i},';'])
    if isfloat(CurField) && ~any(any(CurField))
        C=rmfield(C,MA_MC{i});
    end
end
MA_MC2=fieldnames(C);
for i=1:length(MA_MC2)
    MA_BMW.(MA_MC2{i})=C.(MA_MC2{i});
end
fprintf('\nDone: Model Estimation & Assessment \n')

end