%% ModelFit_BMW
%
% Fit & assess the model designated by the MA_BMW.mat file
% Write results in the original MA_BMW file
% ------------
% Output=ModelFit_BMW(Input)
% Input original MA_BMW.mat (after model definition & configuration),
% write fit results in the original input file.
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Working_Memory_Modeling_Toolbox
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
Output=zeros(Nsubj,1); AIC=zeros(Nsubj,1); BIC=zeros(Nsubj,1);
AICc=zeros(Nsubj,1); DIC=zeros(Nsubj,1); WAIC=zeros(Nsubj,1);
LME=zeros(Nsubj,1);
% SE=zeros(Nsubj,Nss,9);
% Convert model name
ModelName=MN_Convert_BMW(MA_BMW.Model.Model);
fprintf('\nSubject: %s\n','1')
Config_Fit=MA_BMW.Model;
[Param0,Quality]=Mack_Fit(Data{1},Config_Fit,ModelName,MA_BMW.Constraints,MA_BMW.FitOptions);
Ntrial=length(Data{1}.error);
Nparam=length(Param0);
Param=zeros(Nsubj,Nparam);
Param(1,:)=Param0;
Output(1)=Quality.LP;
% SE(1,1:Nss,:)=Quality.SE;
if Nsubj>1
    for subj=2:Nsubj
        fprintf('\nSubject: %s\n',num2str(subj))
        % Fit
        [Param(subj,:),Quality]=Mack_Fit(Data{subj},Config_Fit,ModelName,MA_BMW.Constraints,MA_BMW.FitOptions);
        Output(subj)=Quality.LP;
        %         SE(subj,1:length(unique(Data{subj}.SS)),:)=Quality.SE;
    end
end

%% AIC, AICc, BIC
% Calculate AIC/AICc/BIC by maximum likelihood
for subj=1:Nsubj
    if isfield(MA_BMW.Criteria,'AIC') && MA_BMW.Criteria.AIC==1, AIC(subj)=-2*Output(subj)+2*Nparam; end
    if isfield(MA_BMW.Criteria,'AICc') && MA_BMW.Criteria.AICc==1, AICc(subj)=-2*Output(subj)+2*Nparam+2*Nparam*(Nparam+1)/(Ntrial-Nparam-1); end
    if isfield(MA_BMW.Criteria,'BIC') && MA_BMW.Criteria.BIC==1, BIC(subj)=-2*Output(subj)+Nparam*log(Ntrial); end
end

%% LME
% % Calculate LME based on Riemann sum
% if isfield(MA_BMW.Criteria, 'LME') && MA_BMW.Criteria.LME==1
%     Opt_LME.Nblock=1000;
%     Opt_LME.lb=MA_BMW.Constraints.lb;
%     Opt_LME.ub=MA_BMW.Constraints.ub;
%     Opt_LME.Model=Config_Fit;
%     Opt_LME.Verbosity=1; % Assign 1 to display the result
%     for subj=1:Nsubj
%         LME(subj)=Mack_LME_RS(ModelName, Data{subj}, Opt_LME, Output(subj));
%     end
% end

% Calculate LME based on MCMC

%% DIC, WAIC
% Calculate DIC/WAIC by MCMC
% if MA_BMW.Criteria.DIC==1 || MA_BMW.Criteria.WAIC==1
%     for subj=1:Nsubj
%         Samples=Mack_MCMC(ModelName, Data, Opt_MCMC);
%         if MA_BMW.Criteria.DIC==1
%             DIC(subj)
%         end
%         if MA_BMW.Criteria.WAIC==1
%             WAIC(subj)
%         end
%     end
% end

%% Epilogue
MA_BMW.Param=Param;
if isfield(MA_BMW.Criteria,'Output') && MA_BMW.Criteria.Output==1
    MA_BMW.Output=Output;
end
% MA_BMW.RMSEA=SE;
MA_BMW.AIC=AIC;
MA_BMW.AICc=AICc;
MA_BMW.BIC=BIC;
MA_BMW.DIC=DIC;
MA_BMW.WAIC=WAIC;
MA_BMW.LME=LME;
if ischar(Input)
    Dir=Input(1:end-10);
    cd(Dir)
    save('MA_BMW','MA_BMW');
end
fprintf('\nDone \n')

end