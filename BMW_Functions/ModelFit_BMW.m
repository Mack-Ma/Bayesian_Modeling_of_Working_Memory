%% ModelFit_BMW
%
% Fit & assess the model(s) designated by MA_BMW
% Write results in the original MA_BMW variable
% ------------
% Output=ModelFit_BMW(Input)
%
% Input original MA_BMW (after model definition & configuration),
% write fit results in the original input variable.
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (Mack_ma2018@outlook.com)
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

%% Fitting
Data=MA_BMW.Data;
% Nss=length(unique(Data{1}.SS));
Nsubj=length(Data);
Quality_Fit=cell(Nsubj,1); Output=zeros(Nsubj,1);
LLH=zeros(Nsubj,1); AIC=zeros(Nsubj,1); BIC=zeros(Nsubj,1);
AICc=zeros(Nsubj,1); DIC=zeros(Nsubj,1); DICs=zeros(Nsubj,1); 
WAIC1=zeros(Nsubj,1); WAIC2=zeros(Nsubj,1); LME=zeros(Nsubj,1);
% Convert model name
ModelName=MN_Convert_BMW(MA_BMW.Model.Model);
fprintf('\nSubject: %s\n','1')
Config_Fit=MA_BMW.Model;
[Param0,Quality]=BMW_Fit(Data{1},Config_Fit,ModelName,MA_BMW.Constraints,MA_BMW.FitOptions);
Ntrial=length(Data{1}.error);
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

%% LLH, AIC, AICc, BIC
if any([MA_BMW.Criteria.LLH, MA_BMW.Criteria.AIC, MA_BMW.Criteria.AICc, MA_BMW.Criteria.BIC]==1)
    fprintf('\nNow starting estimate LLH/AIC/AICc/BIC...\n')
    % Calculate LLH/AIC/AICc/BIC based on maximum likelihood
    if isfield(Config_Fit,'Output') && strcmp(Config_Fit.Output,'LLH') % test mode
        for subj=1:Nsubj
            if MA_BMW.Criteria.LLH==1, LLH(subj)=-Output(subj); end
            if MA_BMW.Criteria.AIC==1, AIC(subj)=-2*Output(subj)+2*Nparam; end
            if MA_BMW.Criteria.AICc==1, AICc(subj)=-2*Output(subj)+2*Nparam+2*Nparam*(Nparam+1)/(Ntrial-Nparam-1); end
            if MA_BMW.Criteria.BIC==1, BIC(subj)=-2*Output(subj)+Nparam*log(Ntrial); end
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
                if MA_BMW.Criteria.LLH==1, LLH(subj)=-Output_MLE; end
                if MA_BMW.Criteria.AIC==1, AIC(subj)=-2*Output_MLE+2*Nparam; end
                if MA_BMW.Criteria.AICc==1, AICc(subj)=-2*Output_MLE+2*Nparam+2*Nparam*(Nparam+1)/(Ntrial-Nparam-1); end
                if MA_BMW.Criteria.BIC==1, BIC(subj)=-2*Output_MLE+Nparam*log(Ntrial); end
            end
        end
    end
    fprintf('\nDone!\n')
end

%% LME, DIC, DIC*, WAIC1, WAIC2
% Calculate LME based on MCMC
if any([MA_BMW.Criteria.LME, MA_BMW.Criteria.DIC, MA_BMW.Criteria.DICs, MA_BMW.Criteria.WAIC1, MA_BMW.Criteria.WAIC2]==1)
    fprintf('\nNow start estimate LME/DIC/DIC*/WAIC1/WAIC2...\n')
    if isfield(Quality,'MCMCResult')
        for subj=1:Nsubj
            fprintf('\nSubject: %d',subj)
            Model_MCMC=Config_Fit;
            Model_MCMC.Model=ModelName;
            Model_MCMC.Constraints=MA_BMW.Constraints;
            MCMCcur=Quality_Fit{subj}.MCMCResult;
            if MA_BMW.Criteria.LME==1 
                Method.IC='LME_HarmonicMean';
                LME(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method); 
            end
            if MA_BMW.Criteria.DIC==1 
                Method.IC='DIC';
                DIC(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method); 
            end
            if MA_BMW.Criteria.DICs==1 
                Method.IC='DIC*';
                DICs(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method); 
            end
            if MA_BMW.Criteria.WAIC1==1 
                Method.IC='WAIC1';
                WAIC1(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method); 
            end
            if MA_BMW.Criteria.WAIC2==1 
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
                [MCMCcur,~]=BMW_MCMC(Model_MCMC, Data, Config_MCMC);
                if MA_BMW.Criteria.LME==1
                    Method.IC='LME_HarmonicMean';
                    LME(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if MA_BMW.Criteria.DIC==1
                    Method.IC='DIC';
                    DIC(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if MA_BMW.Criteria.DICs==1
                    Method.IC='DIC*';
                    DICs(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if MA_BMW.Criteria.WAIC1==1
                    Method.IC='WAIC1';
                    WAIC1(subj)=BMW_GetIC_MCMC(MCMCcur,Model_MCMC,Data{subj},Method);
                end
                if MA_BMW.Criteria.WAIC2==1
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
if isfield(MA_BMW.Criteria,'Output') && MA_BMW.Criteria.Output==1
    MA_BMW.Output=Output;
end
MA_BMW.LLH=LLH;
MA_BMW.AIC=AIC;
MA_BMW.AICc=AICc;
MA_BMW.BIC=BIC;
MA_BMW.DIC=DIC;
MA_BMW.DICs=DICs;
MA_BMW.WAIC1=WAIC1;
MA_BMW.WAIC2=WAIC2;
MA_BMW.LME=LME;
if ischar(Input)
    Dir=Input(1:end-10);
    cd(Dir)
    save('MA_BMW','MA_BMW');
end

end