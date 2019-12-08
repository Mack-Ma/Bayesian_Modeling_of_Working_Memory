%% ModelDefinition_BMW
%
% Define models and output MA_BMW.mat file
% -----------------------
% ## Input ##
% Data/Output Dir, Information regarding the model of interest
% ## Output ##
% MA_BMW.mat
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [MA_BMW]=ModelDefinition_BMW(Input)

fprintf('\n---------------------------------------- \n')
fprintf('\nNow start defining model... \n')
MA_BMW.Model=Input.Model;
if ischar(Input.Data)
    RawData=load(Input.Data);
    DataName=fieldnames(RawData);
    eval(['MA_BMW.Data=RawData.',DataName{1},';'])
elseif iscell(Input.Data)
    MA_BMW.Data=Input.Data;
else
    error('Sorry, data file is not valid...')
end
if isfield(Input,'OutputDir')
cd(Input.OutputDir)
save('MA_BMW','MA_BMW');
end
fprintf('\nDone \n')

end
