%% BestModel_Pie
%
% Construct pie charts of best-model ratio 
% -----------------------
% ## Input ##
% 
% ## Output ##
% 
% -----------------------
% Ma, Tianye
% 10/9/2019

function [Fig]=BestModel_Pie_BMW(BM, Info)
Nmodel=length(Info.ModelName);
V_BM=zeros(1,Nmodel);
Nsubj=length(BM);
for i=1:Nmodel
    V_BM(i)=length(find(BM==i));
end
Fig=pie(V_BM/Nsubj);
legend(Info.ModelName)
title('Best Model(',Info.ICName,')')
end