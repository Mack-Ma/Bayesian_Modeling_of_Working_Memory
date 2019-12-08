%% Bayesian Modeling of Working Memory (BMW) Toolbox
%
% Fit and compare models of working memory
% -----------------------
% ## Experimental Paradigm ##
% Continuous Recall/ Change Detection/ Custom
%
% ## Model space ##
% ### Continuous Recall ###
% Standard Mixture(SM)/Item Limit(IL)/Slots-plus-Averaging(SA)/
% Variable Precision(VP)/VP with capacity(VPcap)/
% Equal Precision(EP)/Categorical SA(cSA)/Categorical VP(cVP)/
% Categorical SM(cSM)/Categorical VPcap(cVPcap)
% Derivatives: Swap Rate/Bias/Precision Fluctuation/Bias Fluctuation
% ### Change Detection ###
% Fixed-Capacity(Central-Probe/Full-Display/Single-Probe)/Signal Detection
% Derivatives: Lapse Rate/Ensemble Encoding
%
% ## Criteria for model assessment ##
% Logarithmic Likelihood(LLH)/Akaike Information Criterion(AIC)/
% Modified Akaike Information Criterion(AICc)/Bayesian Information Criterion(BIC)/
% Deviance Information Criterion(DIC)/
% Watanabe-Akaike Information Criterion(WAIC)/
% Log Model Evidence(LME)/
% 2nd-Level Model Frequency/2nd-Level Exceedance Probability
%
% ## Fitting method ##
% Maximum Likelihood(MLE)
%
% ## Optimization algorithm ##
% fmincon/bads/mads/CM-AES/genetic algorithm/
% simulated annealing/bayesopt/fminsearch
% ------------
% Programmed by Ma, Tianye
% Under the guidance of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function BMW(varargin)
clc % Clean the command window
if nargin==0
    close all % Close all the Matlab windows
    clear variables % Clear all the variables in the workspace
    fprintf('Howdy!\n')
    fprintf('\nWelcome you on %s \n',date)
    fprintf('\n### Bayesian Modeling of Working Memory ### \n')
    fprintf('\nVersion: 1.0 (2019) \n')
    fprintf('\nMemory, Attention & Cognition Lab \n')
end
% Welcome screen
if nargin==0
    set(0, 'DefaultFigureMenu','none')
    WelcomeImage=imread('BMW_icon.png');
    imshow(WelcomeImage,'border','tight','initialmagnification','fit');
    ScreenSize=get(0,'screensize');
    LImage=size(WelcomeImage,2)/(size(WelcomeImage,1)+size(WelcomeImage,2))*ScreenSize(3)*3/5;
    HImage=size(WelcomeImage,1)/(size(WelcomeImage,1)+size(WelcomeImage,2))*ScreenSize(3)*3/5;
    set(gcf,'position',[ScreenSize(3)/12, ScreenSize(4)/3, LImage, HImage])
end
% Read BMW path
global bmw
bmw.path=which('BMW');
bmw.path(end-4:end)=[]; % Dir of BMW
% Check optimization toolbox
OptName={'fmincon','fminsearch','BADS','MADS',...
    'Genetic Algorithm','Simulated Annealing','bayesopt','CM-AES'}';
UnavailOpt=[];
switch exist('fmincon','file')
    case 0, UnavailOpt=[UnavailOpt, 1];
end
switch exist('fminsearch','file')
    case 0, UnavailOpt=[UnavailOpt, 2];
end
switch exist('BADS','file')
    case 0, UnavailOpt=[UnavailOpt, 3];
end
switch exist('patternsearch','file')
    case 0, UnavailOpt=[UnavailOpt, 4];
end
switch exist('ga','file')
    case 0, UnavailOpt=[UnavailOpt, 5];
end
switch exist('simulannealbnd','file')
    case 0, UnavailOpt=[UnavailOpt, 6];
end
switch exist('bayesopt','file')
    case 0, UnavailOpt=[UnavailOpt, 7];
end
switch exist('CM-AES','file')
    case 0, UnavailOpt=[UnavailOpt, 8];
end
OptName(UnavailOpt)=[];
bmw.opt=OptName;
% Add path
addpath(bmw.path)
addpath([bmw.path, 'BMW_GUI'])
addpath([bmw.path, 'BMW_Models\Change_Detection'])
addpath([bmw.path, 'BMW_Models\Continuous_Recall'])
addpath([bmw.path, 'BMW_Models\Custom_Models'])
addpath([bmw.path, 'BMW_Functions'])
addpath([bmw.path, 'BMW_Functions\BMW_Visualization'])
set(0, 'DefaultFigureMenu','figure')
if nargin==0
    Panel_Main; % Open main panel
    fprintf('\nType "bmw.opt" to find out the available optimization methods. \n')
    fprintf('\n')
end

end
