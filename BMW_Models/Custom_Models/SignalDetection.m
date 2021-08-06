%% Function: Fit double-threshold dual-process signal detection
%

function Output=SignalDetection(param, Data, Input)
% Specify parameters
dp=param(1); % d-prime (sensitivity)
beta=param(2); % response criterion
Nparam=2;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'Unequal Variance'))
    Nparam=Nparam+1;
    sigma=param(Nparam); % SD for target density
else
    sigma=1;
end
if any(strcmp(Input.Variants,'Oldness Threshold'))
    Nparam=Nparam+1;
    Ro=param(Nparam); % p(always report target)
else
    Ro=0;
end
if any(strcmp(Input.Variants,'Newness Threshold'))
    Nparam=Nparam+1;
    Rn=param(Nparam); % p(always report target)
else
    Rn=0;
end
Ntarget=Data.Ntarget;
Nlure=Data.Nlure;
Nhit=Data.Nhit;
Nfa=Data.Nfa;
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

p_hit=Ro+(1-Ro)*(1-normcdf(beta, dp, sigma)); % hit rate
p_fa=(1-Rn)*(1-normcdf(beta, 0, 1)); % false-alarm rate
LLH_target=-log(binopdf(Nhit, Ntarget, p_hit));
LLH_lure=-log(binopdf(Nfa, Nlure, p_fa));
LLH=LLH_target+LLH_lure; % - log likelihood

% Posterior
LP=-log(Prior)+LLH; % likelihood*prior

if LP==Inf || isnan(LP)
    LP=realmax('double'); % Output should be a real value
end

% PPD
p_LH=[p_hit*ones(Nhit,1); (1-p_hit)*ones(Ntarget-Nhit,1); p_fa*ones(Nlure,1); (1-p_fa)*ones(Nlure-Nfa,1)];

% Decide output
if strcmp(Input.Output,'LP')
    Output=LP;
elseif strcmp(Input.Output,'LLH')
    Output=LLH;
elseif strcmp(Input.Output,'Prior')
    Output=Prior;
elseif strcmp(Input.Output,'LPPD')
    Output=log(p_LH);
elseif strcmp(Input.Output,'All')
    Output.LP=LP;
    Output.LLH=LLH;
    Output.Prior=Prior;
    Output.LPPD=log(p_LH);
end

end