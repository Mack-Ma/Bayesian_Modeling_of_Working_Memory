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
addpath('Your BMW Toolbox Path')
BMW('Silent'); % Add subfolders

% Configuration
Nset=10; % # of datasets
Ntrial=4000; % # of trials per dataset

%% Simulation
Dataset=cell(1,Nset); % Pre-allocation
SS_range=[1 2 4 6]; % Range of set size
% Generate pdf
Data.error_range=-89:1:90; % Axial data
Data.error=0;
Data.SS=reshape(repmat(SS_range,ceil(Ntrial/length(SS_range)),1),[Ntrial,1]);
Input.PDF=1; % Switch on pdf mode
tau=200; % Range of resource allocation variability
kappa=200; % Range of mean precision at set size 1
power=2;
kappa_r=200;
% Generate data
for set=1:Nset
    data_temp.error=zeros(Ntrial,1);
    param=[kappa,tau,power,kappa_r];
    PDF=Variable_Precision(param,Data,Input); % SS-by-error
    freq_temp=0*PDF;
    for i=1:size(PDF,3)
        freq_temp=mnrnd(ceil(Ntrial/length(SS_range)),PDF); % Sampling based on the multinominal distribution
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
    fprintf('\ndataset: %d, kappa: %d, tau: %.2f, power: %.2f, kappa_r: %.2f\n',set,kappa,tau,power,kappa_r) % Progress
    SimData.Param=param;
    SimData.Data=data_temp;
    Dataset{set}=SimData;
    fprintf('\nDataset %d finished.\n',set)
end

%% Circular SD
CSD=zeros(Nset,length(SS_range));
for set=1:Nset
    for ss=1:length(SS_range)
        data_error=Dataset{set}.Data.error;
        data_ss=Dataset{set}.Data.SS;
        data=data_error(data_ss==SS_range(ss));
        CSD(set,ss)=CircSummary_BMW('CircSD',data,180); % Note that the period should be revised accordingly
    end
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
