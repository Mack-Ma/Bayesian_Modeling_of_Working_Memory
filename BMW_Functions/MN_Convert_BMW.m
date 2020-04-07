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
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function [Name]=MN_Convert_BMW(Name0)

% Function names
List_Name={'Item_Limit','Standard_Mixture','Equal_Precision','Slots_plus_Averaging',...
    'Variable_Precision','Variable_Precision_with_Capacity', 'Categorical_Slots_plus_Averaging_BV',...
    'Categorical_Variable_Precision_BV', 'Categorical_Variable_Precision_with_Capacity_BV',...
    'Categorical_Slots_plus_Averaging_WV',...
    'Fixed_Capacity_Central_Probe', 'Fixed_Capacity_Single_Probe', 'Fixed_Capacity_Full_Display',...
    'Signal_Detection'};
% Model names
List_Name0={'Item Limit', 'Standard Mixture', 'Equal Precision', 'Slots-plus-Averaging',...
    'Variable Precision', 'Variable Precision with Capacity', 'Categorical Slots-plus-Averaging (Between-Variant)',...
    'Categorical Variable Precision (Between-Variant)', 'Categorical Variable Precision with Capacity (Between-Variant)',...
    'Categorical Slots-plus-Averaging (Within-Variant)',...
    'Fixed-Capacity (Central Probe)', 'Fixed-Capacity (Single-Probe)', 'Fixed-Capacity (Full-Display)',...
    'Signal Detection'};
Name=List_Name(ismember(List_Name0,Name0));
% For custom models
if isempty(Name)
    Name=Name0;
else
    Name=cell2mat(Name);
end

end