%% Configuration_BMW
%
% Configure model estimation & assessment
% Write information in the original MA_BMW.mat file
% -----------------------
% [MA_BMW]=Configuration_BMW(Input)
% ## Input ##
% Input.InputFile: MA_BMW.mat/directory (after model definition)
% Input.FitOptions, Input.Criteria, Input.Constraints
% ## Output ##
% MA_BMW.mat
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [MA_BMW]=Configuration_BMW(Input)

fprintf('\n---------------------------------------- \n')
fprintf('\nNow start configuration... \n')
if ischar(Input.InputFile)
    load(Input.InputFile);
    DataDir=Input.InputFile;
    Dir=DataDir(1:end-10);
elseif isstruct(Input.InputFile)
    MA_BMW=Input.InputFile;
else
    error('Sorry, the input file is not valid...')
end
% Default set
if isfield(Input.FitOptions,'Default') && Input.FitOptions.Default==1
    Input.FitOptions=Info_BMW('Fit Options'); end
if isfield(Input.Criteria,'Default') && Input.Criteria.Default==1
    Input.Criteria=Info_BMW('Criteria'); end
if isfield(Input.Constraints,'Default') && Input.Constraints.Default==1
    QModel.Model=MA_BMW.Model.Model;
    QModel.Variants=MA_BMW.Model.Variants;
    Input.Constraints=Info_BMW(QModel, MA_BMW.Data);
end
% Output
MA_BMW.FitOptions=ConfigFitOptions_BMW(Input.FitOptions);
MA_BMW.Criteria=Input.Criteria;
MA_BMW.Constraints=Input.Constraints;
if ischar(Input.InputFile)
    cd(Dir)
    save('MA_BMW','MA_BMW');
end
fprintf('\nDone \n')

end