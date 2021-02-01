%% Configure Fit Options
%
% Adjust the input structural array to fit BMW_Fit function.
% ------------
% FitOptions=ConfigFitOptions_BMW(Input)
% Input original FitOptions,
% output converted FitOptions.
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 11/1/2019
%
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/BMW-Ma/Bayesian_Modeling_of_Working_Memory
%

function FitOptions=ConfigFitOptions_BMW(Input)
FitOptions=Input;
if isfield(Input,'Algorithm')
    Algorithm=Input.Algorithm;
else
    Algorithm='GA'; % Default
end
if isfield(Input,'Display')
    Display=Input.Display;
else
    Display='iter'; % Default
end
if isfield(Input,'MaxIter')
    MaxIter=Input.MaxIter;
else
    MaxIter=3000; % Default
end
switch Algorithm
    case 'fmincon: sqp'
        FitOptions.fminconOptions=optimoptions(@fmincon,'Algorithm','sqp');
        % Algorithm: 'active-set'/'interior-point'/'sqp'/'sqp-legacy'/'trust-region-reflective'
        FitOptions.fminconOptions.Display=Display;
        FitOptions.fminconOptions.MaxIter=MaxIter;
        FitOptions.fminconOptions.StepTolerance=1e-6;
    case 'fmincon: interior-point'
        FitOptions.fminconOptions=optimoptions(@fmincon,'Algorithm','interior-point');
        FitOptions.fminconOptions.Display=Display;
        FitOptions.fminconOptions.MaxIter=MaxIter;
        FitOptions.fminconOptions.StepTolerance=1e-6;
    case 'fmincon: active-set'
        FitOptions.fminconOptions=optimoptions(@fmincon,'Algorithm','active-set');
        FitOptions.fminconOptions.Display=Display;
        FitOptions.fminconOptions.MaxIter=MaxIter;
        FitOptions.fminconOptions.StepTolerance=1e-6;
    case 'SA'
        FitOptions.saOptions.Display=Display;
        FitOptions.saOptions.MaxIter=MaxIter;
        FitOptions.saOptions.StepTolerance=1e-6;
    case 'GA'
        FitOptions.gaOptions.Display=Display;
        FitOptions.gaOptions.MaxIter=MaxIter;
        FitOptions.gaOptions.StepTolerance=1e-6;
    case 'MADS'
        FitOptions.madsOptions.Display=Display;
        FitOptions.madsOptions.MaxIter=MaxIter;
        FitOptions.madsOptions.StepTolerance=1e-6;
    case 'BADS'
        FitOptions.badsOptions=bads('defaults');
        FitOptions.badsOptions.Display=Display;
        FitOptions.badsOptions.MaxIter=MaxIter;
        
end
end
