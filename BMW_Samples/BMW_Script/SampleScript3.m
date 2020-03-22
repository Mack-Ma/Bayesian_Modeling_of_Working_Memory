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
% % Load toolbox
% addpath('Your BMW Toolbox Path')
% BMW('Silent'); % Add subfolders
% Designate models based on which to simulate
SimModelSpace={'Item Limit', 'Standard Mixture', 'Slots-plus-Averaging', 'Equal Precision'}; % Full model space
Nsim=length(SimModelSpace);
% Designate models to fit
FitModelSpace={'Item Limit', 'Standard Mixture', 'Slots-plus-Averaging', 'Equal Precision'}; % Full model space
Nfit=length(FitModelSpace);
% Configuration
Nset=10; % # of datasets per model
Ntrial=300; % # of trials per dataset

%% Simulation
fprintf('\nNow start simulation...\n')
Dataset=cell(Nsim,Nset);
SimParamRec=cell(Nsim,Nset);
for sim_ind=1:Nsim
    % Configurations
    SS_range=[1 2 4 6];
    Data.error_range=-89:1:90; % Axial Data
    Data.error=0;
    Data.SS=reshape(repmat(SS_range,ceil(Ntrial/length(SS_range)),1),[Ntrial,1]);
    Input.PDF=1; % Switch on pdf mode (instead of likelihood mode)
    FunctionName=MN_Convert_BMW(SimModelSpace{sim_ind}); % Convert model names to function names
    ModelStruct.Model=SimModelSpace{sim_ind};
    % ModelStruct.Swap=1; % Set model variants
    SimInfo=Info_BMW(ModelStruct,Data); % Get default settings
    SimParam=SimInfo.start;
    % Generate data
    for set=1:Nset
        data_temp.error=zeros(Ntrial,1);
        % Generate pdf
        eval(['PDF=',FunctionName,'(SimParam,Data,Input);']) % SS-by-error
        freq_temp=0*PDF.error;
        for i=1:size(freq_temp,3)
            freq_temp(:,:,i)=mnrnd(ceil(Ntrial/length(SS_range)),PDF.error); % Sampling based on the multinominal distribution
        end
        % Construct data based on the frequency
        data_temp.SS=Data.SS;
        data_temp.error_range=Data.error_range;
        for ss=1:length(SS_range)
            error_temp=zeros(ceil(Ntrial/length(SS_range)),1);
            flag_err=1;
            for err=1:length(Data.error_range)
                if freq_temp(ss,err)~=0
                    error_temp(flag_err:flag_err+freq_temp(ss,err)-1)=Data.error_range(err)*ones(freq_temp(ss,err),1);
                    flag_err=flag_err+freq_temp(ss,err);
                end
            end
            data_temp.error(1+(ss-1)*ceil(Ntrial/length(SS_range)):ss*ceil(Ntrial/length(SS_range)))=error_temp;
        end
        fprintf('\nmodel: %s, dataset: %d\n',SimModelSpace{sim_ind},set) % Progress
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
        Config_MA.Model.Variants.Swap=0; % Swap rate
        Config_MA.Model.Variants.PrecF=0; % Fluctuation of precision
        Config_MA.Model.Variants.BiasF=0; % Fluctuation of bias
        Config_MA.Model.Variants.Bias=0; % Bias
        % Model definition
        MA=ModelDefinition_BMW(Config_MA);
        Config_MA.Criteria.LLH=1; % Assign 1 to calculate log likelihood
        Config_MA.Criteria.AIC=1; % Assign 1 to calculate AIC
        Config_MA.Criteria.AICc=1; % Assign 1 to calculate AICc
        Config_MA.Criteria.BIC=1; % Assign 1 to calculate BIC
        Config_MA.Criteria.DIC=0; % Assign 1 to calculate DIC
        Config_MA.Criteria.WAIC=0; % Assign 1 to calculate WAIC
        Config_MA.Criteria.LME=0; % % Assign 1 to calculate log marginal likelihood
        Config_MA.Constraints.Default=1;
        Config_MA.FitOptions.Default=1;
        %     Config_MA.FitOptions.Algorithm='SA'; % Change optimization algorithm (Default: 'GA')
        %     Config_MA.FitOptions.MaxIter=5000; % Change the max # of iteration (Default: 3000)
        %     Config_MA.FitOptions.Display='off'; % Change the display mode (Default: 'iter')
        Config_MA.InputFile=MA;
        % Configuration
        MA=Configuration_BMW(Config_MA);
        % Estimation
        MA=ModelFit_BMW(MA);
        FitResult{fit_ind,sim_ind}=MA;
    end
end

%% Model Comparison
CompareResult=cell(Nsim,Nset);
Config_MC.Criterion='BIC';
for sim_ind=1:Nsim
    Config_MC.MS=FitResult(:,sim_ind); % Define model space
    CompareResult(sim_ind,:)=ModelComparison_BMW(Config_MC);
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
rmpath('Your Data Path')
