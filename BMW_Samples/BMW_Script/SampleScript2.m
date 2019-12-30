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
Nset=10; % # of datasets per condition
Ntrial=3000; % # of trials per dataset per set size

%% Simulation
Dataset=cell(1,Nset); % Pre-allocation
SS_range=[1 3]; % Range of set size
% Generate pdf
Data.error_range=-89:1:90; % Axial data
Data.error=0;
Data.SS=ones(length(Data.error),1);
Input.PDF=1; % Switch on pdf mode
tau_range=100:20:300; % Range of resource allocation variability
kappa_range=[100 200 300]; % Range of mean precision at set size 1
power=2;
kappa_r=125;
% Generate data
for set=1:Nset
    flag_tau=0;
    SimData=cell(length(tau_range),length(kappa_range));
    for tau=tau_range
        flag_tau=flag_tau+1;
        flag_kappa=0;
        for kappa_bar=kappa_range
            flag_kappa=flag_kappa+1;
            data_temp=zeros(Ntrial,length(SS_range));
            param=[kappa_bar,tau,kappa_r];
            for iSS=1:length(SS_range)
                SS=SS_range(iSS);
                param(1)=param(1)/SS^power;
                Data.SS=SS*Data.SS;
                VPpdf=Variable_Precision(param,Data,Input);
                PDF=[Data.error_range;VPpdf];
                data_temp(:,iSS)=randsrc(Ntrial,1,PDF);
                fprintf('\ndataset: %d, set size: %d, kappa: %.2f, tau: %.2f\n',set,SS,kappa_bar,tau) % Progress
            end
            SimData{flag_tau,flag_kappa}.Param=param;
            SimData{flag_tau,flag_kappa}.Data=data_temp;
        end
    end
    Dataset{set}=SimData;
    fprintf('\nDataset %d finished.\n',set)
end

%% Circular SD
CSD=zeros(Nset,length(tau_range),length(kappa_range),length(SS_range));
for set=1:Nset
    for tau=1:length(tau_range)
        for kappa=1:length(kappa_range)
            for SS=1:length(SS_range)
                data=Dataset{set}{tau,kappa}.Data;
                CSD(set,tau,kappa,SS)=CircSummary_BMW('CircSD',data(:,SS),180); % Note that the period should be revised accordingly
            end
        end
    end
end

%% Relationship between parameters
% CSD-tau-mean precision
for SS=1:length(SS_range)
    figure(SS)
    for kappa=1:length(kappa_range)
        data=CSD(:,:,kappa,SS);
        hold on
        errorbar(tau_range,mean(data),std(data))
    end
end

%% Epilogue
% Save
cd('Your save directory') % cd to the directory where you gonna save the results
save('Your Filename') % Save all variables
% Remove paths
rmpath('Your BMW Toolbox Path')
