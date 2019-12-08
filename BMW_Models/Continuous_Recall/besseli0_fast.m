%% zeroth-order Bessel function
%
% function I0 = besseli0_fast(kappa,scaleflag)
% 
% Returns bessel function for zeroth order, real arguments. Faster than the
% built in matlab function. If scaleflag is set to 1, the output is scaled
% by exp(kappa)
%
% Code copy-pasted from
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/129418,
% which apparently was based on the algorithm given in Numerical Recipes.
%
% This code is part of the paper: 
% - van den Berg, R., Awh, E., & Ma, W. J. (2014). Factorial comparison of working memory models.
% Psychological Review, 121(1), 124.
%

function I0 = besseli0_fast(x,scaleflag)

if ~exist('scaleflag','var')
    scaleflag=0;
end

ax = abs(x);
I0 = zeros(size(x));

% ax<3.75
k = find(ax<3.75);
y=x(k)./3.75;
y=y.^2;
I0(k)=1.0+y.*(3.5156229+y.*(3.0899424+y.*(1.2067492+y.*(0.2659732+y.*(0.360768e-1+y.*0.45813e-2)))));
if scaleflag
    I0(k)=I0(k)./exp(ax(k));
end

% ax>=3.75
k = find(ax>=3.75);
y=3.75./ax(k);
if scaleflag
    I0(k)=1./sqrt(ax(k)).*(0.39894228+y.*(0.1328592e-1+y.*(0.225319e-2+y.*(-0.157565e-2+y.*(0.916281e-2+y.*(-0.2057706e-1+y.*(0.2635537e-1+y.*(-0.1647633e-1+y.*0.392377e-2))))))));
else
    I0(k)=exp(ax(k))./sqrt(ax(k)).*(0.39894228+y.*(0.1328592e-1+y.*(0.225319e-2+y.*(-0.157565e-2+y.*(0.916281e-2+y.*(-0.2057706e-1+y.*(0.2635537e-1+y.*(-0.1647633e-1+y.*0.392377e-2))))))));
end