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
            LLH_BestModel_0=find(LLH_models(subj,:)==max(LLH_models(subj,:)));
            LLH_BestModel(subj)=LLH_BestModel_0(1);
        end
        MC_BMW.LLH.BestModel=LLH_BestModel;
        MC_BMW.LLH.LLH_GLR=LLH_GLR;
        MC_BMW.LLH.delta_LLH=delta_LLH;
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
            AIC_BestModel_0=find(AIC_models(subj,:)==min(AIC_models(subj,:)));
            AIC_BestModel(subj)=AIC_BestModel_0(1);
        end
        MC_BMW.AIC.BestModel=AIC_BestModel;
        MC_BMW.AIC.AIC_GLR=AIC_GLR;
        MC_BMW.AIC.delta_AIC=delta_AIC;
        MC_BMW.AIC.AIC_weight=AIC_weight;
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
            AICc_BestModel_0=find(AICc_models(subj,:)==min(AICc_models(subj,:)));
            AICc_BestModel(subj)=AICc_BestModel_0(1);
        end
        MC_BMW.AICc.BestModel=AICc_BestModel;
        MC_BMW.AICc.AICc_GLR=AICc_GLR;
        MC_BMW.AICc.delta_AICc=delta_AICc;
        MC_BMW.AICc.AICc_weight=AICc_weight;
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
            BIC_BestModel_0=find(BIC_models(subj,:)==min(BIC_models(subj,:)));
            BIC_BestModel(subj)=BIC_BestModel_0(1);
        end
        MC_BMW.BIC.BestModel=BIC_BestModel;
        MC_BMW.BIC.BIC_GBF=BIC_GBF;
        MC_BMW.BIC.delta_BIC=delta_BIC;
        MC_BMW.BIC.BIC_PPr=BIC_PPr;
        % 2nd Level RFX-BMS
        BMC_Results=Mack_BMC(-BIC_models'/2, Opt_BMC);
        MC_BMW.AIC.ModelFreq=BMC_Results.r;
        MC_BMW.AIC.EP=BMC_Results.EP;
    case 'DIC'
        delta_DIC=zeros(Nmodel,Nmodel,Nsubj);
        DIC_PPr=zeros(Nmodel,Nsubj);
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
            DIC_BestModel(subj)=find(DIC_models(subj,:)==min(DIC_models(subj,:)));
        end
        MC_BMW.DIC.BestModel=DIC_BestModel;
        MC_BMW.DIC.delta_DIC=delta_DIC;
        MC_BMW.DIC.DIC_PPr=DIC_PPr;
    case 'DICs'
        delta_DICs=zeros(Nmodel,Nmodel,Nsubj);
        DICs_PPr=zeros(Nmodel,Nsubj);
        DICs_BestModel=zeros(Nsubj,1);
        DICs_models=zeros(Nsubj,Nmodel);
        for subj=1:Nsubj
            for model1=1:Nmodel
                DICs_models(subj,model1)=Q_BMW{model1}.DICs(subj);
                for model2=1:Nmodel
                    delta_DICs(model1,model2,:)=Q_BMW{model2}.DICs-Q_BMW{model1}.DICs;
                end
                delta_DICs_benchmark=delta_DICs(model1,:,subj)-min(delta_DICs(model1,:,subj));
                DICs_PPr(model1,subj)=exp(-delta_DICs_benchmark(model1)/2)/sum(exp(-delta_DICs_benchmark/2));
            end
            DICs_BestModel(subj)=find(DICs_models(subj,:)==min(DICs_models(subj,:)));
        end
        MC_BMW.DICs.BestModel=DICs_BestModel;
        MC_BMW.DICs.delta_DICs=delta_DICs;
        MC_BMW.DICs.DICs_PPr=DICs_PPr;
    case 'WAIC1'
        delta_WAIC1=zeros(Nmodel,Nmodel,Nsubj);
        WAIC1_PPr=zeros(Nmodel,Nsubj);
        WAIC1_BestModel=zeros(Nsubj,1);
        for subj=1:Nsubj
            WAIC1_models=zeros(1,Nmodel);
            for model1=1:Nmodel
                WAIC1_models(model1)=Q_BMW{model1}.WAIC1(subj);
                for model2=1:Nmodel
                    delta_WAIC1(model1,model2,:)=Q_BMW{model2}.WAIC1-Q_BMW{model1}.WAIC1;
                end
                delta_WAIC1_benchmark=delta_WAIC1(model1,:,subj)-min(delta_WAIC1(model1,:,subj));
                WAIC1_PPr(model1,subj)=exp(-delta_WAIC1_benchmark(model1)/2)/sum(exp(-delta_WAIC1_benchmark/2));
            end
            WAIC1_BestModel(subj)=find(WAIC1_models==min(WAIC1_models));
        end
        MC_BMW.WAIC1.BestModel=WAIC1_BestModel;
        MC_BMW.WAIC1.delta_WAIC1=delta_WAIC1;
        MC_BMW.WAIC1.WAIC1_PPr=WAIC1_PPr;
    case 'WAIC2'
        delta_WAIC2=zeros(Nmodel,Nmodel,Nsubj);
        WAIC2_PPr=zeros(Nmodel,Nsubj);
        WAIC2_BestModel=zeros(Nsubj,1);
        for subj=1:Nsubj
            WAIC2_models=zeros(1,Nmodel);
            for model1=1:Nmodel
                WAIC2_models(model1)=Q_BMW{model1}.WAIC2(subj);
                for model2=1:Nmodel
                    delta_WAIC2(model1,model2,:)=Q_BMW{model2}.WAIC2-Q_BMW{model1}.WAIC2;
                end
                delta_WAIC2_benchmark=delta_WAIC2(model1,:,subj)-min(delta_WAIC2(model1,:,subj));
                WAIC2_PPr(model1,subj)=exp(-delta_WAIC2_benchmark(model1)/2)/sum(exp(-delta_WAIC2_benchmark/2));
            end
            WAIC2_BestModel(subj)=find(WAIC2_models==min(WAIC2_models));
        end
        MC_BMW.WAIC2.BestModel=WAIC2_BestModel;
        MC_BMW.WAIC2.delta_WAIC2=delta_WAIC2;
        MC_BMW.WAIC2.WAIC2_PPr=WAIC2_PPr;
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