%% Model Name Conversion
%
% Convert model names into function names
% -----------------------
% FunctionName=MN_Convert_BMW(ModelName)
%
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Name]=MN_Convert_BMW(Name0)

% Function names
List_Name={'Item_Limit','Standard_Mixture','Equal_Precision','Slots_plus_Averaging',...
    'Variable_Precision','Variable_Precision_with_Capacity','Category_Only','Category_Only_with_Capacity',...
    'Fixed_Capacity_Central_Probe', 'Fixed_Capacity_Single_Probe', 'Fixed_Capacity_Full_Display',...
    'Signal_Detection'};
% Model names
List_Name0={'Item Limit', 'Standard Mixture', 'Equal Precision', 'Slots-plus-Averaging',...
    'Variable Precision', 'Variable Precision with Capacity','Category-Only','Category-Only (with Capacity)',...
    'Fixed-Capacity (Central-Probe)', 'Fixed-Capacity (Single-Probe)', 'Fixed-Capacity (Full-Display)',...
    'Signal Detection'};
Name=List_Name(ismember(List_Name0,Name0));
% For custom models
if isempty(Name)
    Name=Name0;
else
    Name=cell2mat(Name);
end

end