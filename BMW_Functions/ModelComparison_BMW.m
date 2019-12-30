%% ModelComparison_BMW
%
% Comparing the models designated by the MA_BMW.mat file
% Write results in the MC_BMW file
% -----------------------
% [MC_BMW]=ModelComparison_BMW(Input)
%
% ## Input ##
% Input.MS: cell array of results of model assessment/ directories of MA_BMW.mat
% Input.Criterion, (Input.MN, Input.OutputDir)
% ## Output ##
% MC_BMW.mat
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [MC_BMW]=ModelComparison_BMW(Input)

fprintf('\n---------------------------------------- \n')
fprintf('\nNow start model comparison... \n')
if ischar(Input.MS)
    DataDir=Input.MS;
    MC_BMW=Input;
    Nmodel=size(DataDir,1);
    load(DataDir(1,:));
    benchmark=MA_BMW;
    Nsubj=size(benchmark.Param,1);
    % Load data
    Q_BMW=cell(1,Nmodel);
    for model=1:Nmodel
        load(DataDir(model,:));
        Q_BMW{model}=MA_BMW;
    end
elseif iscell(Input.MS)
    Nmodel=length(Input.MS);
    Q_BMW=Input.MS;
    Nsubj=size(Q_BMW{1}.Param,1);
else
    error('Sorry, the input file is not valid...')
end

%% Comparison
% Configure BMC
Opt_BMC.Start=ones(Nmodel,1); % Flat prior
Opt_BMC.MaxIter=1500;
Opt_BMC.Stop=1e-6;
Opt_BMC.Verbosity=0;
Opt_BMC.Rec=0;
switch Input.Criterion
    case 'LLH'
        delta_LLH=zeros(Nmodel,Nmodel,Nsubj);
        LLH_GLR=ones(Nmodel,Nmodel);
        LLH_BestModel=zeros(Nsubj,1);
        LLH_models=zeros(Nsubj,Nmodel);
        for subj=1:Nsubj
            for model1=1:Nmodel
                LLH_models(subj,model1)=Q_BMW{model1}.LLH(subj);
                for model2=1:Nmodel
                    delta_LLH(model1,model2,:)=Q_BMW{model2}.LLH-Q_BMW{model1}.LLH;
                end
            end
            LLH_GLR=LLH_GLR.*(exp(LLH_models(subj,:))'*(ones(1,Nmodel)./exp(LLH_models(subj,:))));
            LLH_BestModel(subj)=find(LLH_models(subj,:)==max(LLH_models(subj,:)));
        end
        MC_BMW.LLH.BestModel=LLH_BestModel;
        MC_BMW.LLH.LLH_GLR=LLH_GLR;
        MC_BMW.LLH.delta_LLH=delta_LLH;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(LLH_models', Opt_BMC);
        MC_BMW.LLH.ModelFreq=BMC_Results.r;
        MC_BMW.LLH.EP=BMC_Results.EP;
    case 'AIC'
        delta_AIC=zeros(Nmodel,Nmodel,Nsubj);
        AIC_weight=zeros(Nmodel,Nsubj);
        AIC_GLR=ones(Nmodel,Nmodel);
        AIC_BestModel=zeros(Nsubj,1);
        AIC_models=zeros(Nsubj,Nmodel);
        for subj=1:Nsubj
            for model1=1:Nmodel
                AIC_models(subj,model1)=Q_BMW{model1}.AIC(subj);
                for model2=1:Nmodel
                    delta_AIC(model1,model2,:)=Q_BMW{model2}.AIC-Q_BMW{model1}.AIC;
                end
                delta_AIC_benchmark=delta_AIC(model1,:,subj)-min(delta_AIC(model1,:,subj));
                AIC_weight(model1,subj)=exp(-delta_AIC_benchmark(model1)/2)/sum(exp(-delta_AIC_benchmark/2));
            end
            AIC_GLR=AIC_GLR.*(AIC_weight(:,subj)*(ones(Nmodel,1)./AIC_weight(:,subj))');  % (1,2) = Models1/Model2
            AIC_BestModel(subj)=find(AIC_models(subj,:)==min(AIC_models(subj,:)));
        end
        MC_BMW.AIC.BestModel=AIC_BestModel;
        MC_BMW.AIC.AIC_GLR=AIC_GLR;
        MC_BMW.AIC.delta_AIC=delta_AIC;
        MC_BMW.AIC.AIC_weight=AIC_weight;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(-AIC_models', Opt_BMC);
        MC_BMW.AIC.ModelFreq=BMC_Results.r;
        MC_BMW.AIC.EP=BMC_Results.EP;
    case 'AICc'
        delta_AICc=zeros(Nmodel,Nmodel,Nsubj);
        AICc_weight=zeros(Nmodel,Nsubj);
        AICc_GLR=ones(Nmodel,Nmodel);
        AICc_BestModel=zeros(Nsubj,1);
        AICc_models=zeros(Nsubj,Nmodel);
        for subj=1:Nsubj
            for model1=1:Nmodel
                AICc_models(subj,model1)=Q_BMW{model1}.AICc(subj);
                for model2=1:Nmodel
                    delta_AICc(model1,model2,:)=Q_BMW{model2}.AIC-Q_BMW{model1}.AIC;
                end
                delta_AICc_benchmark=delta_AICc(model1,:,subj)-min(delta_AICc(model1,:,subj));
                AICc_weight(model1,subj)=exp(-delta_AICc_benchmark(model1)/2)/sum(exp(-delta_AICc_benchmark/2));
            end
            AICc_GLR=AICc_GLR.*(AICc_weight(:,subj)*(ones(Nmodel,1)./AICc_weight(:,subj))');
            AICc_BestModel(subj)=find(AICc_models(subj,:)==min(AICc_models(subj,:)));
        end
        MC_BMW.AICc.BestModel=AICc_BestModel;
        MC_BMW.AICc.AICc_GLR=AICc_GLR;
        MC_BMW.AICc.delta_AICc=delta_AICc;
        MC_BMW.AICc.AICc_weight=AICc_weight;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(-AICc_models', Opt_BMC);
        MC_BMW.AICc.ModelFreq=BMC_Results.r;
        MC_BMW.AICc.EP=BMC_Results.EP;
    case 'BIC'
        delta_BIC=zeros(Nmodel,Nmodel,Nsubj);
        BIC_PPr=zeros(Nmodel,Nsubj);
        BIC_GBF=ones(Nmodel,Nmodel);
        BIC_BestModel=zeros(Nsubj,1);
        BIC_models=zeros(Nsubj,Nmodel);
        for subj=1:Nsubj
            for model1=1:Nmodel
                BIC_models(subj,model1)=Q_BMW{model1}.BIC(subj);
                for model2=1:Nmodel
                    delta_BIC(model1,model2,:)=Q_BMW{model2}.BIC-Q_BMW{model1}.BIC;
                end
                delta_BIC_benchmark=delta_BIC(model1,:,subj)-min(delta_BIC(model1,:,subj));
                % Consider a flat prior
                BIC_PPr(model1,subj)=exp(-delta_BIC_benchmark(model1)/2)/sum(exp(-delta_BIC_benchmark/2));
            end
            BIC_GBF=BIC_GBF.*(BIC_PPr(:,subj)*(ones(Nmodel,1)./BIC_PPr(:,subj))');
            BIC_BestModel(subj)=find(BIC_models(subj,:)==min(BIC_models(subj,:)));
        end
        MC_BMW.BIC.BestModel=BIC_BestModel;
        MC_BMW.BIC.BIC_GBF=BIC_GBF;
        MC_BMW.BIC.delta_BIC=delta_BIC;
        MC_BMW.BIC.BIC_PPr=BIC_PPr;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(-BIC_models', Opt_BMC);
        MC_BMW.BIC.ModelFreq=BMC_Results.r;
        MC_BMW.BIC.EP=BMC_Results.EP;
    case 'DIC'
        delta_DIC=zeros(Nmodel,Nmodel,Nsubj);
        DIC_PPr=zeros(Nmodel,Nsubj);
        DIC_GBF=ones(Nmodel,Nmodel);
        DIC_BestModel=zeros(Nsubj,1);
        DIC_models=zeros(Nsubj,Nmodel);
        for subj=1:Nsubj
            for model1=1:Nmodel
                DIC_models(subj,model1)=Q_BMW{model1}.DIC(subj);
                for model2=1:Nmodel
                    delta_DIC(model1,model2,:)=Q_BMW{model2}.DIC-Q_BMW{model1}.DIC;
                end
                delta_DIC_benchmark=delta_DIC(model1,:,subj)-min(delta_DIC(model1,:,subj));
                DIC_PPr(model1,subj)=exp(-delta_DIC_benchmark(model1)/2)/sum(exp(-delta_DIC_benchmark/2));
            end
            DIC_GBF=DIC_GBF.*(DIC_PPr(:,subj)*(ones(Nmodel,1)./DIC_PPr(:,subj))');
            DIC_BestModel(subj)=find(DIC_models(subj,:)==min(DIC_models(subj,:)));
        end
        MC_BMW.DIC.BestModel=DIC_BestModel;
        MC_BMW.DIC.DIC_GBF=DIC_GBF;
        MC_BMW.DIC.delta_DIC=delta_DIC;
        MC_BMW.DIC.DIC_PPr=DIC_PPr;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(DIC_models', Opt_BMC);
        MC_BMW.DIC.ModelFreq=BMC_Results.r;
        MC_BMW.DIC.EP=BMC_Results.EP;
    case 'WAIC'
        delta_WAIC=zeros(Nmodel,Nmodel,Nsubj);
        WAIC_PPr=zeros(Nmodel,Nsubj);
        WAIC_GBF=ones(Nmodel,Nmodel);
        WAIC_BestModel=zeros(Nsubj,1);
        for subj=1:Nsubj
            WAIC_models=zeros(1,Nmodel);
            for model1=1:Nmodel
                WAIC_models(model1)=Q_BMW{model1}.WAIC(subj);
                for model2=1:Nmodel
                    delta_WAIC(model1,model2,:)=Q_BMW{model2}.WAIC-Q_BMW{model1}.WAIC;
                end
                delta_WAIC_benchmark=delta_WAIC(model1,:,subj)-min(delta_WAIC(model1,:,subj));
                WAIC_PPr(model1,subj)=exp(-delta_WAIC_benchmark(model1)/2)/sum(exp(-delta_WAIC_benchmark/2));
            end
            WAIC_GBF=WAIC_GBF.*(WAIC_PPr(:,subj)*(ones(Nmodel,1)./WAIC_PPr(:,subj))');
            WAIC_BestModel(subj)=find(WAIC_models==min(WAIC_models));
        end
        MC_BMW.WAIC.BestModel=WAIC_BestModel;
        MC_BMW.WAIC.WAIC_GBF=WAIC_GBF;
        MC_BMW.WAIC.delta_WAIC=delta_WAIC;
        MC_BMW.WAIC.WAIC_PPr=WAIC_PPr;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(MC_BMW.WAIC.WAIC_PPr, Opt_BMC);
        MC_BMW.WAIC.ModelFreq=BMC_Results.r;
        MC_BMW.WAIC.EP=BMC_Results.EP;
    case 'LME'
        delta_LME=zeros(Nmodel,Nmodel,Nsubj);
        LME_Group=zeros(Nmodel,Nsubj);
        LME_GBF=ones(Nmodel,Nmodel);
        LME_BestModel=zeros(Nsubj,1);
        for subj=1:Nsubj
            LME_models=zeros(subj,Nmodel);
            for model1=1:Nmodel
                LME_models(subj,model1)=Q_BMW{model1}.LME(subj);
                for model2=1:Nmodel
                    delta_LME(model1,model2,:)=Q_BMW{model2}.LME-Q_BMW{model1}.LME;
                end
            end
            LME_GBF=LME_GBF.*(LME_models(subj,:)'*(ones(1,Nmodel)./LME_models(subj,:)));
            LME_BestModel(subj)=find(LME_models(subj,:)==max(LME_models(subj,:)));
            LME_Group(:,subj)=LME_models(subj,:)';
        end
        MC_BMW.LME.BestModel=LME_BestModel;
        MC_BMW.LME.LME_GBF=LME_GBF;
        MC_BMW.LME.delta_LME=delta_LME;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(LME_Group, Opt_BMC);
        MC_BMW.LME.ModelFreq=BMC_Results.r;
        MC_BMW.LME.EP=BMC_Results.EP;
end

%% Epilogue
if isfield(Input,'MN')
    MC_BMW.Names=Input.MN;
end
if isfield(Input,'OutputDir')
    cd(Input.OutputDir)
    save('MC_BMW','MC_BMW');
end
fprintf('\nDone \n')

end