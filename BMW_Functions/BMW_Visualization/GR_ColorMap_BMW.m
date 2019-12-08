%% GR_ColorMap
%
% Construct colormaps of log group bayesian factor/ log IC ratio 
% -----------------------
% ## Input ##
% 
% ## Output ##
% 
% -----------------------
% Ma, Tianye
% 10/9/2019

function [Fig]=GR_ColorMap_BMW(GR, Info)

Nmodel=size(GR,1);

Fig=figure;

% Create axes
axes1 = axes('Parent',Fig);
hold(axes1,'on');

% Create image
image(GR,'Parent',axes1,'CDataMapping','scaled');
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'Layer','top','XTick',1:Nmodel,'XTickLabel',Info.ModelName,...
    'YTick',1:Nmodel,'YTickLabel',Info.ModelName);
% Create colorbar
cb=colorbar('peer',axes1);
colorTitleHandle = get(cb,'Title');
set(colorTitleHandle ,'String',Info.ICName);

end