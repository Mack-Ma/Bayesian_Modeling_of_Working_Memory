%% Calculate LME (Riemann Sum)
%
% Calculate log model evidence by Riemann sum
% -----------------------
% ## Input ##
% model, Data, config, CenteredValue
% ## Output ##
% LME
% ## Reference ##
% - Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., & Friston, K. J. (2009). 
% "Bayesian model selection for group studies". NeuroImage, 46(4), 1004-1017.
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/27/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox:
% https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function LME=Mack_LME_RS(model, Data, config, CenteredValue)
Nblock=config.Nblock;
lb=config.lb;
ub=config.ub;
verb=config.Verbosity;
Config_Model=config.Model;
Nparam=length(lb);
step=zeros(1,Nparam);
for i=1:Nparam
    step(i)=(ub(i)-lb(i))/Nblock;
end
% Uniform prior
Prior=1/(Nblock).^Nparam;
Area_Step=1;
for i=1:Nparam
    Area_Step0=step(i);
    Area_Step=Area_Step*Area_Step0;
end
% Calculate LME by Riemann sum
ME=0;
param=lb+step; % Initial value
seq=1:Nparam;
count=0;
if verb==1
    fprintf('\nNow start calculating LME by Riemann sum... \n')
end
while 1
    count=count+1;
    if count>=ceil((Nblock*Nparam)/10) && verb==1
        fprintf('|||||') % Progress bar
        count=1;
    end
    eval(['LLH=-',model,'(param, Data, Config_Model);'])
    ME=ME+exp(LLH+log(Prior)+log(Area_Step)-CenteredValue); % Due to computational limit
    i=randsample(seq,1);
    % Random walk
    param(i)=param(i)+step(i);
    if param(i)>=ub(i) && i<=length(seq)
        seq(i)=[];
    end
    if isempty(seq)
        LME=log(ME)+CenteredValue;
        break
    end
end
if verb==1
    fprintf(' Done! \nLME=%.2d\n',LME)
end
end
