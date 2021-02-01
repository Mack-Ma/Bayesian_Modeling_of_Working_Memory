%% CircSummary_BMW
%
% Calculate circular summary statistics
% ------------
% Output=CircSummary_BMW(Stat, Data, Period)
%
% - Stat
% A string that specifies the output, potential choices are,
%       'Mean', mean resultant vector direction in degree
%       'NormSD', normal standard deviation based on circular difference
%       'Moment', 1st trigonometric moment about the zero direction
%       'AngVar', circular variance
%       'AngDev', mean angular deviation
%       'CircSD', circular standard deviation
%       'CircDev', circular mean deviation
%       'CosDist', mean cosine distance
%       'Kurtosis', circular kurtosis
%       'Skewness', circular skewness
%
% - Data
% An N-element vector that contains angular data in degree
%
% - Period
% A constant in degree that specifies the circular period (e.g. 180),
% i.e. -Period/2 < Data <= Period/2,
% default as 360
%
% - Bin
% An integer that specifies the number of bins which determines the correction factor,
% default as Inf
%
% ------------
% Programmed by Ma, Tianye
% Mentored by Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 12/9/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function Output=CircSummary_BMW(Stat, Data, Period, config)

if nargin<4 || (nargin==4 && ~isfield(config,'Bin'))
    c1=1; % 1st correction factor
    c2=1; % 2nd correction factor
    if nargin==2
        Period=360;
    end
else
    Bin=config.Bin;
    UnitBin=Period/Bin;
    c1=UnitBin/2/sind(UnitBin/2);
    c2=UnitBin/sind(UnitBin);
end

if nargin<4 || (nargin==4 && ~isfield(config,'weight'))
    config.weight=ones(size(Data));
end
weight=config.weight/sum(config.weight);

Order=360/Period;
A=mean(weight.*cosd(Order*Data));
B=mean(weight.*sind(Order*Data));
Moment=A+1i*B; % 1st Trigonometric moment
R=c1*sqrt(A^2+B^2); % Mean resultant vector length
MeanAng=atan2d(B, A)/Order; % Mean resultant vector direction
V=1-R; % Angular variance
AD=sqrt(2*(1-R)); % Angular deviation
CSD=sqrt(-2*log(R)); % Circular standard deviation

if strcmp(Stat,'Mean')
    Output=MeanAng;
elseif strcmp(Stat, 'NormSD')
    Output=sqrt(mean((CircDist_BMW('Ang',Data,MeanAng)).^2)); % Normal SD
elseif strcmp(Stat,'Moment')
    Output=Moment;
elseif strcmp(Stat,'AngVar')
    Output=V*180/pi;
elseif strcmp(Stat, 'AngDev')
    Output=AD*180/pi;
elseif strcmp(Stat,'CircSD')
    Output=CSD*180/pi;
elseif strcmp(Stat,'CircDev')
    if size(Data,2)==1 % When Data is a column vector
        M1=repmat(Data, 1, length(Data));
        M2=repmat(Data', length(Data), 1);
    elseif size(Data,1)==1 % When Data is a row vector
        M1=repmat(Data, length(Data), 1);
        M2=repmat(Data', 1, length(Data));
    end
    CMD=mean(mean(CircDist_BMW('Ang',M1,M2,Period)));
    Output=CMD;
elseif strcmp(Stat, 'CosDist')
    CosD=1-R^2; % Mean cosine distance
    Output=CosD;
elseif strcmp(Stat, 'Dispersion')
    A2=mean(cosd(2*Order*Data));
    B2=mean(sind(2*Order*Data));
    R2=c2*sqrt(A2^2+B2^2); % Length of second mean resultant vector
    Disp=(1-R2)/2*R^2; % Circular dispersion
    Output=Disp;
elseif strcmp(Stat, 'Skewness')
    A2=mean(cosd(2*Order*Data));
    B2=mean(sind(2*Order*Data));
    R2=c2*sqrt(A2^2+B2^2); % Length of second mean resultant vector
    Ang2=atan2(B2,A2)/Order; % Direction of second mean resultant vector
    SinM2=R2*sind(Ang2-2*MeanAng); % 2nd central sine moment
    S=SinM2/(1-R)^1.5; % Skewness
    Output=S;
elseif strcmp(Stat, 'Kurtosis')
    A2=mean(cosd(2*Order*Data));
    B2=mean(sind(2*Order*Data));
    R2=c2*sqrt(A2^2+B2^2); % Length of second mean resultant vector
    Ang2=atan2(B2,A2)/Order; % Direction of second mean resultant vector
    CosM2=R2*cosd(Ang2-2*MeanAng); % 2nd central sine moment
    K=(CosM2-R^4)/(1-R)^2; % Kurtosis
    Output=K;
end

end