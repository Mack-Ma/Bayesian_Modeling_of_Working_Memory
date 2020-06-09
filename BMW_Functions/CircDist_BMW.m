%% CircDist_BMW
%
% Calculate circular distance
% ------------
% Output=CircDist_BMW(Stat, Data1, Data2, Period)
%
% - Stat
% A string that specifies the output, potential choices are,
%       'Diff', angular difference (with directionality: Data1-Data2)
%       'Cos', cosine distance
%
% - Data1
% An N-element vector that contains angular data in degree
%
% - Data2
% An N-element vector or a constant in degree
%
% - Period
% A constant in degree that specifies the circular period (e.g. 180),
% default as 360
%
% ------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 12/9/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function Output=CircDist_BMW(Stat, Data1, Data2, Period)
if nargin==3
    Period=360;
end
if length(Data2)==1 % When Data2 is a constant
    Data2=repmat(Data2,size(Data1,1),size(Data1,2));
end
N=length(Data1);
if N~=length(Data2)
    error('Data1 and Data2 should have the same length...')
end
if strcmp(Stat, 'Ang')
    Output=Period/2-abs(Period/2-abs(Data1-Data2));
elseif strcmp(Stat, 'Diff')
    Output=mod(Data1-Data2+Period/2,Period)-Period/2;
elseif strcmp(Stat, 'Cos')
    Output=1-cosd(Data1-Data2);
end
end