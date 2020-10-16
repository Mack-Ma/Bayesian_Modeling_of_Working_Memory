%% ModelComparison_BMW
%
% Comparing models based on the designated model criterion
% -----------------------
% MC_BMW=ModelComparison_BMW(Input)
%
% For detailed information, plz refer to the manual (type BMW('Manual') in the command line)
% For the complete example, plz refer to "SampleScript1" in the toolbox
% ## Input ##
% - Input 
%   cell array, with each element as the output of ModelFit_BMW for each
%   model.
%
% ## Output ##
% - MC_BMW
%   struct, with each subfield comprises the results of the comparison
%   based on each model criterion
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [MC_BMW]=ModelComparison_BMW(Input, family)

fprintf('\n---------------------------------------- \n')
fprintf('\nNow start model comparison... \n')
if ischar(Input)
    DataDir=Input;
    RawMC=load(DataDir);
    MC_Field=fieldnames(RawMC);
    if length(MC_Field)>1
        error('Sorry, the input file is not valid...')
    end
    eval(['MC=RawMC.',MC_Field{1},';'])
    if length(size(MC))==2 && any(size(MC)==1)
        Nmodel=length(MC);
        Nrepeat=1;
    else
        Nrepeat=size(MC,1);
        Nmodel=size(MC,2);
    end
    Nsubj=size(MC{1}.Param,1);
    Q_BMW0=MC;
elseif iscell(Input)
    Q_BMW0=Input;
     if length(size(Input))==2 && size(Input,1)==1
        Nmodel=length(Input);
        Nrepeat=1;
     elseif length(size(Input))==2 && size(Input,1)~=1
        Nrepeat=size(Input,1);
        Nmodel=size(Input,2);
    end
    Nsubj=size(Q_BMW0{1}.Param,1);
else
    error('Sorry, the input file is not valid...')
end
if nargin==1
    family_bool=0;
else
    family_bool=1;
    family_range=unique(family(family~=0));
    Nfamily=length(family_range);
end

%% Comparison
% Configure BMS
MC_BMW=cell(Nrepeat,1);
for r=1:Nrepeat
    Q_BMW=Q_BMW0(r,:);
    if Q_BMW{1}.FitOptions.RFXBMS==1
        Opt_BMC.Start=ones(Nmodel,1); % Flat prior
        Opt_BMC.MaxIter=5000;
        Opt_BMC.Stop=1e-6;
        Opt_BMC.Verbosity=0;
        Opt_BMC.Rec=0;
        if family_bool==1
            Opt_BMC_fam=Opt_BMC;
            Opt_BMC_fam.Start=ones(Nfamily,1);
        end
    end
    for i=1:length(Q_BMW{1}.Criteria)
        switch Q_BMW{1}.Criteria{i}
            case 'LLH'
                delta_LLH=zeros(Nmodel,Nmodel,Nsubj);
                LLH_GLR=ones(Nmodel,Nmodel);
                LLH_weight=zeros(Nmodel,Nsubj);
                LLH_BestModel=zeros(Nsubj,1);
                LLH_models=zeros(Nsubj,Nmodel);
                for subj=1:Nsubj
                    for model1=1:Nmodel
                        LLH_models(subj,model1)=Q_BMW{model1}.LLH(subj);
                        sum_LH=0;
                        for model2=1:Nmodel
                            delta_LLH(model1,model2,:)=Q_BMW{model2}.LLH-Q_BMW{model1}.LLH;
                            sum_LH=sum_LH+exp(Q_BMW{model2}.LLH(subj));
                        end
                        delta_LLH_benchmark=delta_LLH(model1,:,subj)-min(delta_LLH(model1,:,subj));
                        LLH_weight(model1,subj)=exp(-delta_LLH_benchmark(model1))/sum(exp(-delta_LLH_benchmark));
                    end
                    LLH_GLR=LLH_GLR.*(exp(LLH_models(subj,:))'*(ones(1,Nmodel)./exp(LLH_models(subj,:))));
                    LLH_BestModel_0=find(LLH_models(subj,:)==max(LLH_models(subj,:)));
                    LLH_BestModel(subj)=LLH_BestModel_0(1);
                end
                MC_BMW{r}.LLH.BestModel=LLH_BestModel;
                MC_BMW{r}.LLH.LLH_GLR=LLH_GLR;
                MC_BMW{r}.LLH.delta_LLH=delta_LLH;
                MC_BMW{r}.LLH.LLH_weight=LLH_weight';
                allcumsum_LLH=cumsum(LLH_weight');
                MC_BMW{r}.LLH.Group_PP=allcumsum_LLH(end,:)/sum(allcumsum_LLH(end,:));
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
                MC_BMW{r}.AIC.BestModel=AIC_BestModel;
                MC_BMW{r}.AIC.AIC_GLR=AIC_GLR;
                MC_BMW{r}.AIC.delta_AIC=delta_AIC;
                MC_BMW{r}.AIC.AIC_weight=AIC_weight';
                allcumsum_AIC=cumsum(AIC_weight');
                MC_BMW{r}.AIC.Group_PP=allcumsum_AIC(end,:)/sum(allcumsum_AIC(end,:));
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
                MC_BMW{r}.AICc.BestModel=AICc_BestModel;
                MC_BMW{r}.AICc.AICc_GLR=AICc_GLR;
                MC_BMW{r}.AICc.delta_AICc=delta_AICc;
                MC_BMW{r}.AICc.AICc_weight=AICc_weight';
                allcumsum_AICc=cumsum(AICc_weight');
                MC_BMW{r}.AICc.Group_PP=allcumsum_AICc(end,:)/sum(allcumsum_AICc(end,:));
            case 'BIC'
                if family_bool==1
                    delta_BIC_Family=zeros(Nfamily,Nfamily,Nsubj);
                    BIC_weight_Family=zeros(Nfamily,Nsubj);
                    BIC_GBF_Family=ones(Nfamily,Nfamily);
                    BIC_BestFamily=zeros(Nfamily,1);
                    BIC_Family=zeros(Nsubj,Nfamily);
                end
                delta_BIC=zeros(Nmodel,Nmodel,Nsubj);
                BIC_weight=zeros(Nmodel,Nsubj);
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
                        BIC_weight(model1,subj)=exp(-delta_BIC_benchmark(model1)/2)/sum(exp(-delta_BIC_benchmark/2));
                    end
                    BIC_GBF=BIC_GBF.*(BIC_weight(:,subj)*(ones(Nmodel,1)./BIC_weight(:,subj))');
                    BIC_BestModel_0=find(BIC_models(subj,:)==min(BIC_models(subj,:)));
                    BIC_BestModel(subj)=BIC_BestModel_0(1);
                    if family_bool==1
                        for family1=1:Nfamily
                            fam_ind=family==family_range(family1);
                            BIC_Family(subj,family1)=-2*log(mean(exp(-(BIC_models(subj,fam_ind)-max(BIC_models(subj,fam_ind)))/2)))+max(BIC_models(subj,fam_ind));
                            for family2=1:Nfamily
                                delta_BIC_Family(family1,family2,:)=BIC_Family(:,family2)-BIC_Family(:,family1);
                            end
                            delta_BIC_benchmark_fam=delta_BIC_Family(family1,:,subj)-min(delta_BIC_Family(family1,:,subj));
                            BIC_weight_Family(family1,subj)=exp(-delta_BIC_benchmark_fam(family1)/2)/sum(exp(-delta_BIC_benchmark_fam/2));
                        end
                        BIC_GBF_Family=BIC_GBF_Family.*(BIC_weight_Family(:,subj)*(ones(Nfamily,1)./BIC_weight_Family(:,subj))');
                        BIC_BestFamily_0=find(BIC_Family(subj,:)==min(BIC_Family(subj,:)));
                        BIC_BestFamily(subj)=BIC_BestFamily_0(1);
                    end
                end
                MC_BMW{r}.BIC.BestModel=BIC_BestModel;
                MC_BMW{r}.BIC.BIC_GBF=BIC_GBF;
                MC_BMW{r}.BIC.delta_BIC=delta_BIC;
                MC_BMW{r}.BIC.BIC_weight=BIC_weight';
                allcumsum_BIC=cumsum(BIC_weight');
                MC_BMW{r}.BIC.Group_PP=allcumsum_BIC(end,:)/sum(allcumsum_BIC(end,:));
                if family_bool==1
                    MC_BMW{r}.BIC.BestFamily=BIC_BestFamily;
                    MC_BMW{r}.BIC.BIC_GBF_Family=BIC_GBF_Family;
                    MC_BMW{r}.BIC.delta_BIC_Family=delta_BIC_Family;
                    MC_BMW{r}.BIC.BIC_weight_Family=BIC_weight_Family';
                    allcumsum_BIC_fam=cumsum(BIC_weight_Family');
                    MC_BMW{r}.BIC.Group_PP_Family=allcumsum_BIC_fam(end,:)/sum(allcumsum_BIC_fam(end,:));
                end
                % 2nd Level RFX-BMS
                if Q_BMW{1}.FitOptions.RFXBMS==1
                    BMC_Results=BMW_BMS(-BIC_models'/2, Opt_BMC);
                    MC_BMW{r}.BIC.ModelFreq=BMC_Results.r;
                    MC_BMW{r}.BIC.EP=BMC_Results.EP;
                    if family_bool==1
                        BMC_Results_fam=BMW_BMS(-BIC_Family'/2, Opt_BMC_fam);
                        MC_BMW{r}.BIC.FamilyFreq=BMC_Results_fam.r;
                        MC_BMW{r}.BIC.EP_Family=BMC_Results_fam.EP;
                    end
                end
            case 'DIC1'
                delta_DIC1=zeros(Nmodel,Nmodel,Nsubj);
                DIC1_weight=zeros(Nmodel,Nsubj);
                DIC1_BestModel=zeros(Nsubj,1);
                DIC1_models=zeros(Nsubj,Nmodel);
                for subj=1:Nsubj
                    for model1=1:Nmodel
                        DIC1_models(subj,model1)=Q_BMW{model1}.DIC1(subj);
                        for model2=1:Nmodel
                            delta_DIC1(model1,model2,:)=Q_BMW{model2}.DIC1-Q_BMW{model1}.DIC1;
                        end
                        delta_DIC1_benchmark=delta_DIC1(model1,:,subj)-min(delta_DIC1(model1,:,subj));
                        DIC1_weight(model1,subj)=exp(-delta_DIC1_benchmark(model1)/2)/sum(exp(-delta_DIC1_benchmark/2));
                    end
                    DIC1_BestModel(subj)=find(DIC1_models(subj,:)==min(DIC1_models(subj,:)));
                end
                MC_BMW{r}.DIC1.BestModel=DIC1_BestModel;
                MC_BMW{r}.DIC1.delta_DIC1=delta_DIC1;
                MC_BMW{r}.DIC1.DIC1_weight=DIC1_weight';
                allcumsum_DIC1=cumsum(DIC1_weight');
                MC_BMW{r}.DIC1.Group_PP=allcumsum_DIC1(end,:)/sum(allcumsum_DIC1(end,:));
            case 'DIC2'
                delta_DIC2=zeros(Nmodel,Nmodel,Nsubj);
                DIC2_weight=zeros(Nmodel,Nsubj);
                DIC2_BestModel=zeros(Nsubj,1);
                DIC2_models=zeros(Nsubj,Nmodel);
                for subj=1:Nsubj
                    for model1=1:Nmodel
                        DIC2_models(subj,model1)=Q_BMW{model1}.DIC2(subj);
                        for model2=1:Nmodel
                            delta_DIC2(model1,model2,:)=Q_BMW{model2}.DIC2-Q_BMW{model1}.DIC2;
                        end
                        delta_DIC2_benchmark=delta_DIC2(model1,:,subj)-min(delta_DIC2(model1,:,subj));
                        DIC2_weight(model1,subj)=exp(-delta_DIC2_benchmark(model1)/2)/sum(exp(-delta_DIC2_benchmark/2));
                    end
                    DIC2_BestModel(subj)=find(DIC2_models(subj,:)==min(DIC2_models(subj,:)));
                end
                MC_BMW{r}.DIC2.BestModel=DIC2_BestModel;
                MC_BMW{r}.DIC2.delta_DIC2=delta_DIC2;
                MC_BMW{r}.DIC2.DIC2_weight=DIC2_weight';
                allcumsum_DIC2=cumsum(DIC2_weight');
                MC_BMW{r}.DIC2.Group_PP=allcumsum_DIC2(end,:)/sum(allcumsum_DIC2(end,:));
            case 'DICs'
                delta_DICs=zeros(Nmodel,Nmodel,Nsubj);
                DICs_weight=zeros(Nmodel,Nsubj);
                DICs_BestModel=zeros(Nsubj,1);
                DICs_models=zeros(Nsubj,Nmodel);
                for subj=1:Nsubj
                    for model1=1:Nmodel
                        DICs_models(subj,model1)=Q_BMW{model1}.DICs(subj);
                        for model2=1:Nmodel
                            delta_DICs(model1,model2,:)=Q_BMW{model2}.DICs-Q_BMW{model1}.DICs;
                        end
                        delta_DICs_benchmark=delta_DICs(model1,:,subj)-min(delta_DICs(model1,:,subj));
                        DICs_weight(model1,subj)=exp(-delta_DICs_benchmark(model1)/2)/sum(exp(-delta_DICs_benchmark/2));
                    end
                    DICs_BestModel(subj)=find(DICs_models(subj,:)==min(DICs_models(subj,:)));
                end
                MC_BMW{r}.DICs.BestModel=DICs_BestModel;
                MC_BMW{r}.DICs.delta_DICs=delta_DICs;
                MC_BMW{r}.DICs.DICs_weight=DICs_weight';
                allcumsum_DICs=cumsum(DICs_weight');
                MC_BMW{r}.DICs.Group_PP=allcumsum_DICs(end,:)/sum(allcumsum_DICs(end,:));
            case 'WAIC1'
                delta_WAIC1=zeros(Nmodel,Nmodel,Nsubj);
                WAIC1_weight=zeros(Nmodel,Nsubj);
                WAIC1_BestModel=zeros(Nsubj,1);
                for subj=1:Nsubj
                    WAIC1_models=zeros(1,Nmodel);
                    for model1=1:Nmodel
                        WAIC1_models(model1)=Q_BMW{model1}.WAIC1(subj);
                        for model2=1:Nmodel
                            delta_WAIC1(model1,model2,:)=Q_BMW{model2}.WAIC1-Q_BMW{model1}.WAIC1;
                        end
                        delta_WAIC1_benchmark=delta_WAIC1(model1,:,subj)-min(delta_WAIC1(model1,:,subj));
                        WAIC1_weight(model1,subj)=exp(-delta_WAIC1_benchmark(model1)/2)/sum(exp(-delta_WAIC1_benchmark/2));
                    end
                    WAIC1_BestModel(subj)=find(WAIC1_models==min(WAIC1_models));
                end
                MC_BMW{r}.WAIC1.BestModel=WAIC1_BestModel;
                MC_BMW{r}.WAIC1.delta_WAIC1=delta_WAIC1;
                MC_BMW{r}.WAIC1.WAIC1_weight=WAIC1_weight';
                allcumsum_WAIC1=cumsum(WAIC1_weight');
                MC_BMW{r}.WAIC1.Group_PP=allcumsum_WAIC1(end,:)/sum(allcumsum_WAIC1(end,:));
            case 'WAIC2'
                delta_WAIC2=zeros(Nmodel,Nmodel,Nsubj);
                WAIC2_weight=zeros(Nmodel,Nsubj);
                WAIC2_BestModel=zeros(Nsubj,1);
                for subj=1:Nsubj
                    WAIC2_models=zeros(1,Nmodel);
                    for model1=1:Nmodel
                        WAIC2_models(model1)=Q_BMW{model1}.WAIC2(subj);
                        for model2=1:Nmodel
                            delta_WAIC2(model1,model2,:)=Q_BMW{model2}.WAIC2-Q_BMW{model1}.WAIC2;
                        end
                        delta_WAIC2_benchmark=delta_WAIC2(model1,:,subj)-min(delta_WAIC2(model1,:,subj));
                        WAIC2_weight(model1,subj)=exp(-delta_WAIC2_benchmark(model1)/2)/sum(exp(-delta_WAIC2_benchmark/2));
                    end
                    WAIC2_BestModel(subj)=find(WAIC2_models==min(WAIC2_models));
                end
                MC_BMW{r}.WAIC2.BestModel=WAIC2_BestModel;
                MC_BMW{r}.WAIC2.delta_WAIC2=delta_WAIC2;
                MC_BMW{r}.WAIC2.WAIC2_weight=WAIC2_weight';
                allcumsum_WAIC2=cumsum(WAIC2_weight');
                MC_BMW{r}.WAIC2.Group_PP=allcumsum_WAIC2(end,:)/sum(allcumsum_WAIC2(end,:));
            case 'LME_BS'
                if family_bool==1
                    delta_LME_BS_Family=zeros(Nfamily,Nfamily,Nsubj);
                    LME_BS_weight_Family=zeros(Nfamily,Nsubj);
                    LME_BS_GBF_Family=ones(Nfamily,Nfamily);
                    LME_BS_BestFamily=zeros(Nfamily,1);
                    LME_BS_Family=zeros(Nsubj,Nfamily);
                end
                delta_LME_BS=zeros(Nmodel,Nmodel,Nsubj);
                LME_BS_Group=zeros(Nmodel,Nsubj);
                LME_BS_GBF=ones(Nmodel,Nmodel);
                LME_BS_BestModel=zeros(Nsubj,1);
                LME_BS_weight=zeros(Nmodel,Nsubj);
                for subj=1:Nsubj
                    LME_BS_models=zeros(subj,Nmodel);
                    for model1=1:Nmodel
                        LME_BS_models(subj,model1)=-Q_BMW{model1}.LME_BS(subj);
                        for model2=1:Nmodel
                            delta_LME_BS(model1,model2,:)=Q_BMW{model2}.LME_BS-Q_BMW{model1}.LME_BS;
                        end
                        delta_LME_BS_benchmark=delta_LME_BS(model1,:,subj)-min(delta_LME_BS(model1,:,subj));
                        LME_BS_weight(model1,subj)=exp(-delta_LME_BS_benchmark(model1))/sum(exp(-delta_LME_BS_benchmark));
                    end
                    LME_BS_GBF=LME_BS_GBF.*(LME_BS_models(subj,:)'*(ones(1,Nmodel)./LME_BS_models(subj,:)));
                    LME_BS_BestModel(subj)=find(LME_BS_models(subj,:)==max(LME_BS_models(subj,:)));
                    LME_BS_Group(:,subj)=LME_BS_models(subj,:)';
                    if family_bool==1
                        for family1=1:Nfamily
                            fam_ind=family==family_range(family1);
                            LME_BS_Family(subj,family1)=-log(mean(exp(-(LME_BS_models(subj,fam_ind)-max(LME_BS_models(subj,fam_ind))))))+max(LME_BS_models(subj,fam_ind));
                            for family2=1:Nfamily
                                delta_LME_BS_Family(family1,family2,:)=LME_BS_Family(:,family2)-LME_BS_Family(:,family1);
                            end
                            delta_LME_BS_benchmark_fam=delta_LME_BS_Family(family1,:,subj)-min(delta_LME_BS_Family(family1,:,subj));
                            LME_BS_weight_Family(family1,subj)=exp(-delta_LME_BS_benchmark_fam(family1))/sum(exp(-delta_LME_BS_benchmark_fam));
                        end
                        LME_BS_GBF_Family=LME_BS_GBF_Family.*(LME_BS_weight_Family(:,subj)*(ones(Nfamily,1)./LME_BS_weight_Family(:,subj))');
                        LME_BS_BestFamily_0=find(LME_BS_Family(subj,:)==min(LME_BS_Family(subj,:)));
                        LME_BS_BestFamily(subj)=LME_BS_BestFamily_0(1);
                    end
                end
                MC_BMW{r}.LME_BS.BestModel=LME_BS_BestModel;
                MC_BMW{r}.LME_BS.LME_BS_GBF=LME_BS_GBF;
                MC_BMW{r}.LME_BS.delta_LME_BS=delta_LME_BS;
                MC_BMW{r}.LME_BS.LME_BS_weight=LME_BS_weight';
                allcumsum_LME_BS=cumsum(LME_BS_weight');
                MC_BMW{r}.LME_BS.Group_PP=allcumsum_LME_BS(end,:)/sum(allcumsum_LME_BS(end,:));
                if family_bool==1
                    MC_BMW{r}.LME_BS.BestFamily=LME_BS_BestFamily;
                    MC_BMW{r}.LME_BS.LME_BS_GBF_Family=LME_BS_GBF_Family;
                    MC_BMW{r}.LME_BS.delta_LME_BS_Family=delta_LME_BS_Family;
                    MC_BMW{r}.LME_BS.LME_BS_weight_Family=LME_BS_weight_Family';
                    allcumsum_LME_BS_fam=cumsum(LME_BS_weight_Family');
                    MC_BMW{r}.LME_BS.Group_PP_Family=allcumsum_LME_BS_fam(end,:)/sum(allcumsum_LME_BS_fam(end,:));
                end
                % 2nd Level RFX-BMS
                if Q_BMW{1}.FitOptions.RFXBMS==1
                    BMC_Results=BMW_BMS(LME_BS_Group, Opt_BMC);
                    MC_BMW{r}.LME_BS.ModelFreq=BMC_Results.r;
                    MC_BMW{r}.LME_BS.EP=BMC_Results.EP;
                    if family_bool==1
                        BMC_Results_fam=BMW_BMS(LME_BS_Family', Opt_BMC_fam);
                        MC_BMW{r}.LME_BS.FamilyFreq=BMC_Results_fam.r;
                        MC_BMW{r}.LME_BS.EP_Family=BMC_Results_fam.EP;
                    end
                end
            case 'LME_GHM'
                if family_bool==1
                    delta_LME_GHM_Family=zeros(Nfamily,Nfamily,Nsubj);
                    LME_GHM_weight_Family=zeros(Nfamily,Nsubj);
                    LME_GHM_GBF_Family=ones(Nfamily,Nfamily);
                    LME_GHM_BestFamily=zeros(Nfamily,1);
                    LME_GHM_Family=zeros(Nsubj,Nfamily);
                end
                delta_LME_GHM=zeros(Nmodel,Nmodel,Nsubj);
                LME_GHM_Group=zeros(Nmodel,Nsubj);
                LME_GHM_GBF=ones(Nmodel,Nmodel);
                LME_GHM_BestModel=zeros(Nsubj,1);
                LME_GHM_weight=zeros(Nmodel,Nsubj);
                for subj=1:Nsubj
                    LME_GHM_models=zeros(subj,Nmodel);
                    for model1=1:Nmodel
                        LME_GHM_models(subj,model1)=-Q_BMW{model1}.LME_GHM(subj);
                        for model2=1:Nmodel
                            delta_LME_GHM(model1,model2,:)=Q_BMW{model2}.LME_GHM-Q_BMW{model1}.LME_GHM;
                        end
                        delta_LME_GHM_benchmark=delta_LME_GHM(model1,:,subj)-min(delta_LME_GHM(model1,:,subj));
                        LME_GHM_weight(model1,subj)=exp(-delta_LME_GHM_benchmark(model1))/sum(exp(-delta_LME_GHM_benchmark));
                    end
                    LME_GHM_GBF=LME_GHM_GBF.*(LME_GHM_models(subj,:)'*(ones(1,Nmodel)./LME_GHM_models(subj,:)));
                    LME_GHM_BMlist=find(LME_GHM_models(subj,:)==max(LME_GHM_models(subj,:)));
                    LME_GHM_BestModel(subj)=LME_GHM_BMlist(randsample(1:length(LME_GHM_BMlist),1));
                    LME_GHM_Group(:,subj)=LME_GHM_models(subj,:)';
                    if family_bool==1
                        for family1=1:Nfamily
                            fam_ind=family==family_range(family1);
                            LME_GHM_Family(subj,family1)=-log(mean(exp(-(LME_GHM_models(subj,fam_ind)-max(LME_GHM_models(subj,fam_ind))))))+max(LME_GHM_models(subj,fam_ind));
                            for family2=1:Nfamily
                                delta_LME_GHM_Family(family1,family2,:)=LME_GHM_Family(:,family2)-LME_GHM_Family(:,family1);
                            end
                            delta_LME_GHM_benchmark_fam=delta_LME_GHM_Family(family1,:,subj)-max(delta_LME_GHM_Family(family1,:,subj));
                            LME_GHM_weight_Family(family1,subj)=exp(delta_LME_GHM_benchmark_fam(family1))/sum(exp(delta_LME_GHM_benchmark_fam));
                        end
                        LME_GHM_GBF_Family=LME_GHM_GBF_Family.*(LME_GHM_weight_Family(:,subj)*(ones(Nfamily,1)./LME_GHM_weight_Family(:,subj))');
                        LME_GHM_BestFamily_0=find(LME_GHM_Family(subj,:)==min(LME_GHM_Family(subj,:)));
                        LME_GHM_BestFamily(subj)=LME_GHM_BestFamily_0(1);
                    end
                end
                MC_BMW{r}.LME_GHM.BestModel=LME_GHM_BestModel;
                MC_BMW{r}.LME_GHM.LME_GHM_GBF=LME_GHM_GBF;
                MC_BMW{r}.LME_GHM.delta_LME_GHM=delta_LME_GHM;
                MC_BMW{r}.LME_GHM.LME_GHM_weight=LME_GHM_weight';
                allcumsum_LME_GHM=cumsum(LME_GHM_weight');
                MC_BMW{r}.LME_GHM.Group_PP=allcumsum_LME_GHM(end,:)/sum(allcumsum_LME_GHM(end,:));
                if family_bool==1
                    MC_BMW{r}.LME_GHM.BestFamily=LME_GHM_BestFamily;
                    MC_BMW{r}.LME_GHM.LME_GHM_GBF_Family=LME_GHM_GBF_Family;
                    MC_BMW{r}.LME_GHM.delta_LME_GHM_Family=delta_LME_GHM_Family;
                    MC_BMW{r}.LME_GHM.LME_GHM_weight_Family=LME_GHM_weight_Family';
                    allcumsum_LME_GHM_fam=cumsum(LME_GHM_weight_Family');
                    MC_BMW{r}.LME_GHM.Group_PP_Family=allcumsum_LME_GHM_fam(end,:)/sum(allcumsum_LME_GHM_fam(end,:));
                end
                % 2nd Level RFX-BMS
                if Q_BMW{1}.FitOptions.RFXBMS==1
                    BMC_Results=BMW_BMS(LME_GHM_Group, Opt_BMC);
                    MC_BMW{r}.LME_GHM.ModelFreq=BMC_Results.r;
                    MC_BMW{r}.LME_GHM.EP=BMC_Results.EP;
                    if family_bool==1
                        BMC_Results_fam=BMW_BMS(LME_GHM_Family', Opt_BMC_fam);
                        MC_BMW{r}.LME_GHM.FamilyFreq=BMC_Results_fam.r;
                        MC_BMW{r}.LME_GHM.EP_Family=BMC_Results_fam.EP;
                    end
                end
        end
    end
end

%% Epilogue
fprintf('\nDone: Model Comparison & Selection \n')

end