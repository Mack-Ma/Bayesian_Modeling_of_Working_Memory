%% Convert parameters to fit MCMC
%
% Do transform between the totally real number field and a given real number subfield
%

function NewParam=MCMCConvert_BMW(Param, UB, LB, Method)
% resize
if size(Param,1)~=1
    UB=repmat(UB,[size(Param,1),1]);
    LB=repmat(LB,[size(Param,1),1]);
end
Range=UB-LB;
% set default method
if nargin<3
    Method='InverseFisher';
end
switch Method
    case 'Probit' % [lb, ub]=>(-Inf, +Inf)
        % first project from the given interval to [-1,1]
        NewParam0=(Param-LB)*2./Range-1;
        % do Probit Transform
        NewParam=erfinv(NewParam0);
    case 'InverseProbit' % (-Inf, +Inf)=>[lb, ub]
        % do Transform
        NewParam0=erf(Param);
        % project parameters from [-1,1] to the given interval
        NewParam=(NewParam0+1).*Range/2+LB;
    case 'Logit' % [lb,ub]=>(-Inf, +Inf)
        % first do preliminary projection
        NewParam0=(Param-LB)./Range;
        % do Logit Transform
        NewParam=log(NewParam0)-log(1-NewParam0);
    case 'InverseLogit' % (-Inf, +Inf)=>[lb,ub]
        % check bound
        Param(Param<-700)=-700;
        % do Inverse-Logit Transform
        NewParam0=1./(1+exp(-Param));
        % further scaling
        NewParam=NewParam0*Range-LB;
    case 'Fisher' % [lb, ub] => (-Inf, +Inf) 
        % first project from the given interval to [-1,1]
        NewParam0=(Param-LB)*2./Range-1;
        % do Fisher Transform
        NewParam=0.5*log((1+NewParam0)./(1-NewParam0));
    case 'InverseFisher' % (-Inf, +Inf) => [lb, ub]
        % check range due to the computational limit
        Param(Param>350)=350;
        % do Inverse Fisher Transform
        NewParam0=(exp(2*Param)-1)./(exp(2*Param)+1);
        % project parameters from [-1,1] to the given interval
        NewParam=(NewParam0+1).*Range/2+LB;
    case 'Gaussian' % (-Inf, +Inf) => [lb, ub]
        % do Transform
        NewParam0=normpdf(Param,0,1);
        % project from (0,1) to the given interval
        NewParam=NewParam0*Range+LB;
end
end
