%% Sample Script 3: Model Recovery Analysis
%
% Test the criteria by comparing models using the simulated data
% ------------
% Sample script for BMW toolbox
% For detailed instruction of this script, please refer to the manual (type BMW('Manual')).
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
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
SimModelSpace={'Slots-plus-Averaging','Variable Precision','Categorical Slots-plus-Averaging (Between-Variant)','Categorical Variable Precision (Between-Variant)',...
    'Categorical Slots-plus-Averaging (Within-Variant)','Categorical Variable Precision (Within-Variant)'}; 
Nsim=length(SimModelSpace);
% Designate models to fit
FitModelSpace={'Slots-plus-Averaging','Variable Precision','Categorical Slots-plus-Averaging (Between-Variant)','Categorical Variable Precision (Between-Variant)',...
    'Categorical Slots-plus-Averaging (Within-Variant)','Categorical Variable Precision (Within-Variant)'}; 
Nfit=length(FitModelSpace);
% Configuration
Nset=80; % # of datasets per model
Ntrial=3000; % # of trials per dataset
% Set parameters
load('SimParam_BMW.mat')
Param_Sim=cell(1,Nsim);
for i=1:Nsim
    par=SimParam_BMW{i};
    Param_Sim{i}=repmat(par,ceil(Nset/size(par,1)),1);
end

%% Simulation
fprintf('\nNow start simulation...\n')
Dataset=cell(Nsim,Nset);
SimParamRec=cell(Nsim,Nset);
for sim_ind=1:Nsim
    % Configurations
    SS_range=[1 2 4 6];
    Data.sample_range=1:180;
    Data.response_range=1:180;
    Data.category_range=[90 180];
    Data.error_range=-90:89;
    for set=1:Nset
        % Generate samples
        samples_ss1=randsample(Data.sample_range,30*length(Data.response_range))';
        % Get categories
        cat_diff=abs(CircDist_BMW('Diff',repmat(Data.category_range,[length(samples_ss1),1]),repmat(samples_ss1,[1,length(Data.category_range)])));
        [~,cat_ind]=min(cat_diff,[],2);
        categories_ss1=Data.category_range(cat_ind)';
        samples=repmat(samples_ss1,length(SS_range),1);
        categories=repmat(categories_ss1,length(SS_range),1);
        responses=repmat(repmat(Data.response_range',30,1),length(SS_range),1);
        Data.SS=reshape(repmat(SS_range,30*length(Data.response_range),1),[length(SS_range)*30*length(Data.response_range),1]);
        Input.Output='LPPD'; % use log pointwise predictive density
        FunctionName=MN_Convert_BMW(SimModelSpace{sim_ind}); % Convert model names to function names
        ModelStruct.Model=SimModelSpace{sim_ind};
        % ModelStruct.Variants={'Swap'}; % Set model variants
        SimInfo=Info_BMW(ModelStruct,Data); % Get default settings
        Param_Model=Param_Sim{sim_ind};
        data_temp.error=zeros(Ntrial,1);
        SimParam=Param_Model(set,:);
        % Generate categorical errors
        Data.error_c=randsample(error_c_range,length(errors),true);
        % Generate pdf
        eval(['PDF=',FunctionName,'(SimParam,Data,Input);']) % SS-by-error
        PDF=exp(PDF);
        PDF=reshape(PDF,[length(Data.error_range),length(SS_range)])';
        PDF=PDF./repmat(sum(PDF,2),1,size(PDF,2)); % to ensure that PDF is a valid prob. density function
        freq_temp=mnrnd(ceil(Ntrial/length(SS_range)),PDF); % Sampling based on the multinominal distribution
        % Construct data based on the frequency
        data_temp.SS=reshape(repmat(SS_range,ceil(Ntrial/length(SS_range)),1),[Ntrial,1]);
        data_temp.error_range=Data.error_range;
        data_temp.error_c=Data.error_c;
        for ss=1:length(SS_range)
            error_temp=zeros(ceil(Ntrial/length(SS_range)),1);
            flag_err=1;
            for err=1:length(Data.error_range)
                if freq_temp(ss,err)~=0
                    error_temp(flag_err:flag_err+freq_temp(ss,err)-1)=Data.error_range(err)*ones(freq_temp(ss,err),1);
                    flag_err=flag_err+freq_temp(ss,err);
                end
            end
            error_temp;
            data_temp.error(1+(ss-1)*ceil(Ntrial/length(SS_range)):ss*ceil(Ntrial/length(SS_range)))=error_temp;
        end
        fprintf('\nModel: %s, Dataset: %d\n',SimModelSpace{sim_ind},set) % Progress
        Dataset{sim_ind,set}=data_temp;
        SimParamRec{sim_ind,set}=SimParam;
    end
end

%% Fit
fprintf('\nNow start model fitting...\n')
FitResult=cell(Nfit,Nsim);
for sim_ind=1:Nsim
    for fit_ind=1:Nfit
        fprintf('\n%s\n',FitModelSpace{fit_ind})
        Config_MA.Data=Dataset(sim_ind,:);
        % Specify model
        Config_MA.Model.Model=FitModelSpace{fit_ind};
        % Set Variants (Optional)
%         Config_MA.Model.Variants={'Bias'};
        % Model definition
        Config_MA.Criteria={'DIC1','DIC2','DIC*','WAIC1','WAIC2','LME_GHM','LME_BS'};
        %     Config_MA.FitOptions.Algorithm='SA'; % Change optimization algorithm (Default: 'GA')
        % Configuration
        MA=Configuration_BMW(Config_MA);
        % Estimation
        MA=ModelFit_BMW(MA);
        FitResult{fit_ind,sim_ind}=MA;
    end
end

%% Model Comparison
CompareResult=cell(Nsim,Nfit);
for sim_ind=1:Nsim
    CompareResult(sim_ind,:)=ModelComparison_BMW(FitResult(:,sim_ind));
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
rmpath('Your Data Path')
