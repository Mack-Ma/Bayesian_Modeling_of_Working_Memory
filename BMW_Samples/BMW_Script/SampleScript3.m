%% Sample Script 3: Model Recovery Analysis
%
% Test the criteria by comparing models using the simulated data
% ------------
% Sample script for BMW toolbox
% For detailed instruction of this script, please refer to the manual (type BMW('Manual')).
% ------------
% Programmed by Ma, Tianye
% Memory, Attention & Cognition (MAC) Lab,
% Sun Yat-Sen University
% 1/25/2020
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

%% Prologue
close all
clear variables
% load toolbox
BMW('AddPath');
% Designate models based on which to simulate
SimModelSpace={'Slots-plus-Averaging','Variable Precision'};
SimVariantSpace={{'ContinuousK'},{'ContinuousK','Category (Between-Item)'}, {'ContinuousK','Category (Within-Item)'}};
Nsim_models=length(SimModelSpace);
Nsim_variants=length(SimVariantSpace);
Nsim=Nsim_models*Nsim_variants;
% Designate models to fit
FitModelSpace={'Slots-plus-Averaging','Variable Precision'};
FitVariantSpace={{'ContinuousK'},{'ContinuousK', 'Category (Between-Item)'}, {'ContinuousK','Category (Within-Item)'}};
Nfit_models=length(FitModelSpace);
Nfit_variants=length(FitVariantSpace);
Nfit=Nfit_models*Nfit_variants;
% Configuration
Nset=32; % # of datasets per model
Ntrial=1000; % # of trials per dataset
% Set parameters
load('SimParam_BMW.mat')
Param_Sim=cell(1,Nsim);
for i=1:Nfit
    par=SimParam{i};
    Param_Sim{i}=repmat(par,ceil(Nset/size(par,1)),1);
end

%% Simulation
fprintf('\nNow start simulation...\n')
Dataset=cell(Nsim,Nset);
SimParamRec=cell(Nsim,Nset);
sim_ind=0;
for var_ind=1:Nsim_variants
    for model_ind=1:Nsim_models
        sim_ind=sim_ind+1;
        % Configurations
        SS_range=[1 2 4 6];
        Data.sample_range=1:180;
        Data.response_range=1:180;
        Data.category_range=[90 180];
        Data.error_range=-90:89;
        for set=1:Nset
            % Get pdf
            PDF=zeros(length(SS_range),length(Data.response_range),length(Data.sample_range));
            for s=1:length(Data.sample_range)
                Data.sample=Data.sample_range(s)*ones(length(Data.response_range)*length(SS_range),1);
                % Find the category
                cat_diff0=abs(CircDist_BMW('Diff',repmat(Data.category_range,[length(Data.sample),1]),repmat(Data.sample,[1,length(Data.category_range)])));
                [~,cat_ind0]=min(cat_diff0,[],2);
                Data.category=Data.category_range(cat_ind0)';
                Data.response=repmat(Data.response_range',length(SS_range),1);
                Data.SS=reshape(repmat(SS_range,length(Data.response_range),1),length(Data.response_range)*length(SS_range),1);
                Input.Output='LPPD'; % use log pointwise predictive density
                FunctionName=MN_Convert_BMW(SimModelSpace{model_ind}); % Convert model names to function names
                ModelFunc=str2func(FunctionName);
                Input.Variants=SimVariantSpace{var_ind}; % Set model variants
                Param_Model=Param_Sim{sim_ind};
                Sim_Param=Param_Model(set,:);
                PDF0=exp(ModelFunc(Sim_Param,Data,Input));
                PDF0=reshape(PDF0,[length(Data.error_range),length(SS_range)])';
                PDF0=PDF0./repmat(sum(PDF0,2),1,size(PDF0,2)); % to ensure that PDF is a valid prob. density function
                PDF(:,:,s)=PDF0;
            end
            % Simulation
            % Generate samples
            samples=randsample(Data.sample_range,Ntrial,true)';
            % Get categories
            cat_diff=abs(CircDist_BMW('Diff',repmat(Data.category_range,[length(samples),1]),repmat(samples,[1,length(Data.category_range)])));
            [~,cat_ind]=min(cat_diff,[],2);
            categories=Data.category_range(cat_ind)';
            % Construct data based on the frequency
            data_temp.sample=samples;
            data_temp.category=categories;
            data_temp.SS=reshape(repmat(SS_range,ceil(Ntrial/length(SS_range)),1),[Ntrial,1]);
            data_temp.error_range=Data.error_range;
            data_temp.sample_range=Data.sample_range;
            data_temp.response_range=Data.response_range;
            data_temp.category_range=Data.category_range;
            data_temp.response=zeros(Ntrial,1);
            for trial=1:ceil(Ntrial)
                s=samples(trial);
                ss=SS_range==data_temp.SS(trial);
                response_temp=zeros(ceil(Ntrial/length(SS_range)),1);
                    PDF_now=PDF(ss,:,Data.sample_range==s);
                    freq_temp=mnrnd(1,PDF_now); % Sampling based on the multinominal distribution
                    data_temp.response(trial)=Data.response_range(freq_temp==1);
            end
            if isempty(cell2mat(SimVariantSpace{var_ind}))
                Variants_Display='None';
            else
                Variants_Now=SimVariantSpace{var_ind};
                Variants_Display=cell(1,length(Variants_Now));
                for i=1:length(Variants_Now)
                    Variants_Display{i}=[Variants_Now{i},' '];
                end
                Variants_Display=cell2mat(Variants_Display);
            end
            fprintf('\nModel: %s, Variants: %s, Dataset: %d\n',SimModelSpace{model_ind},Variants_Display,set) % Progress
            Dataset{sim_ind,set}=data_temp;
            SimParamRec{sim_ind,set}=Sim_Param;
        end
    end
end

%% Fit
fprintf('\nNow start model fitting...\n')
FitResult=cell(Nfit,Nsim);
for sim_ind=1:Nsim
    fit_ind=0;
    Config_MA.Data=Dataset(sim_ind,:);
    for var_ind=1:Nfit_variants
        for model_ind=1:Nfit_models
            fit_ind=fit_ind+1;
            if isempty(cell2mat(FitVariantSpace{var_ind}))
                Variants_Display='None';
            else
                Variants_Now=FitVariantSpace{var_ind};
                Variants_Display=cell(1,length(Variants_Now));
                for i=1:length(Variants_Now)
                    Variants_Display{i}=[Variants_Now{i},' '];
                end
                Variants_Display=cell2mat(Variants_Display);
            end
            fprintf('\nModel: %s, Variants: %s\n',FitModelSpace{model_ind},Variants_Display) % Progress
            % Specify model
            Config_MA.Model.Model=FitModelSpace{model_ind};
            % Specify Variants
            Config_MA.Model.Variants=FitVariantSpace{var_ind};
            % Specify the criteria of interest
            Config_MA.Criteria={'DIC1','DIC2','DIC*','WAIC1','WAIC2'};
            %     Config_MA.FitOptions.Algorithm='GA'; % Change the optimization algorithm (Default: 'DE-MCMC')
            % Configuration
            MA=Configuration_BMW(Config_MA);
            % Estimation
            MA=ModelFit_BMW(MA);
            FitResult{sim_ind,fit_ind}=MA;
        end
    end
end

%% Model Comparison
CompareResult=cell(Nsim,1);
for sim_ind=1:Nsim
    CompareResult(sim_ind)=ModelComparison_BMW(FitResult(sim_ind,:));
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
rmpath('Your Data Path')
