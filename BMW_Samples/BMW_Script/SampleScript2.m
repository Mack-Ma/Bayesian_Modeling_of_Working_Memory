%% Sample Script 2: Simulation
%
% Simulate continuous recall data based on Variable Precision model
% ------------
% Sample script for BMW toolbox
% For detailed instruction of this script, please refer to the manual (type BMW('Manual')).
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 12/5/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

%% Prologue
close all
clear variables
SimModel='Variable Precision'; % The model for simulation
% Load toolbox
addpath('Your Toolbox Path')
BMW('AddPath'); % Add subfolders
% Configuration
Nset=10; % # of datasets (per group of parameters)
Ntrial=4000; % # of trials per dataset

%% Simulation
SS_range=[1 2 4 6]; % Range of set size
% Generate pdf
Data.error_range=-89:1:90; % Axial data
Data.error=repmat(Data.error_range,1,length(SS_range));
Data.SS=reshape(repmat(SS_range,length(Data.error_range),1),[length(SS_range)*length(Data.error_range),1]);
Input.Output='LPPD'; % Switch on pdf mode
tau_range=20:10:200; % resource allocation variability
kappa=100; % mean precision at set size 1
power=2;
kappa_r=100;
Dataset=cell(Nset,length(tau_range)); % Pre-allocation
% Generate data
for set=1:Nset
    for t=1:length(tau_range)
        tau=tau_range(t);
        data_temp.error=zeros(Ntrial,1);
        param=[kappa,tau,power,kappa_r];
        PDF=Variable_Precision(param,Data,Input); % SS-by-error
        PDF=exp(PDF);
        PDF=reshape(PDF,[length(SS_range),length(Data.error_range)]);
        PDF=PDF./repmat(sum(PDF,2),1,size(PDF,2)); % to ensure that PDF is a valid prob. density function
        freq_temp=mnrnd(ceil(Ntrial/length(SS_range)),PDF); % Sampling based on the multinominal distribution
        % Construct data based on the frequency
        data_temp.SS=reshape(repmat(SS_range,ceil(Ntrial/length(SS_range)),1),[Ntrial,1]);
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
        fprintf('\ndataset: %d, kappa: %d, tau: %.2f, power: %.2f, kappa_r: %.2f\n',set,kappa,tau,power,kappa_r) % Progress
        SimData.Param=param;
        SimData.Data=data_temp;
        Dataset{set,t}=SimData;
        fprintf('\nDataset %d finished.\n',set)
    end
end

%% Circular SD
CSD=zeros(Nset,length(tau_range),length(SS_range));
for set=1:Nset
    for t=1:length(tau_range)
        for ss=1:length(SS_range)
            data_error=Dataset{set,t}.Data.error;
            data_ss=Dataset{set,t}.Data.SS;
            data=data_error(data_ss==SS_range(ss));
            CSD(set,t,ss)=CircSummary_BMW('CircSD',data,180); % Note that the period should be revised accordingly, help CircSummary_BMW for details
        end
    end
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your Toolbox Path')
