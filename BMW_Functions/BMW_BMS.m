%% Second-Level Bayesian Model Comparison
%
% 2nd-level hierachical model describing the probability of "choosing" a certain model
% based on the model evidence for each model and each subject.
% -----------------------
% Output=BMW_BMC(LME,Config)
%
% ## Input ##
% - LME
% List of log model evidence (rows for each model & cols for each subject)
% - Config
% Options that configure the comparsion process
% ## Output ##
% - output
% Estimated parameters and estimation records of the 2nd-level model
% ## Reference ##
% - Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., & Friston, K. J. (2009). 
% "Bayesian model selection for group studies". NeuroImage, 46(4), 1004-1017.
% -----------------------
% Programmed by Ma, Tianye
% Under the instruction of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% 10/2/2019
%
% Bug reports or any other feedbacks please contact M.T. (BMW_ma2018@outlook.com)
% BMW toolbox: https://github.com/BMW-Ma/Bayesian_Modeling_of_Working_Memory
%

function [output]=BMW_BMS(LME,Config)

if isfield(Config, 'Start'), a0=Config.Start; 
else a0=ones(size(LME,1),1); end % Start values of the parameters of the dirichlet prior
MaxIter=Config.MaxIter; % Maximum times of iteration
Stop=Config.Stop; % Stop criterion
Rec=Config.Rec; % Assign 1 to record the result of each iteration
Verbosity=Config.Verbosity; % Assign 1 to display iterations
Nmodel=size(LME,1); % Total # of models

% Start iteration (variational bayes approximation)
fprintf('\nNow start 2nd-level bayesian model comparison... \n')
a=a0; % Specify start value
FlagIter=0; % # of iteration
RecIter=zeros(1,Nmodel);
while 1
    FlagIter=FlagIter+1; % # of the current iteration
    log_u=LME+repmat(psi(a),1,size(LME,2))-psi(sum(a));
    u=exp(log_u-repmat(max(log_u,[],1),size(LME,1),1)); % Exponentiate
    b=sum(u./repmat(sum(u,1),size(u,1),1),2);
    a_prev=a;
    a=a0+b; % Add "data counts" to the "prior counts"
    if Rec==1
        RecIter(FlagIter,:)=a;
    end
    if norm(a-a_prev)<Stop % Converge
        fprintf('\n\n----------\n')
        fprintf('\n## Break iteration: Converge at %s times of iteration ##\n', num2str(FlagIter))
        fprintf('\n----------\n\n')
        break;
    end
    if FlagIter==MaxIter % Fail convergence after max # of iteration
        fprintf('\n\n----------\n')
        fprintf('\n## Break iteration: Reach max number of iteration (%s) ##\n', num2str(MaxIter))
        fprintf('\n----------\n\n')
        break;
    end
    if Verbosity==1
        fprintf('\n\n----------\n')
        fprintf('\nIter: %s diff: %d \n',num2str(FlagIter),norm(a-a_prev));
        fprintf('\nParameters: \n')
        fprintf('\n %.2d ',a)
    end
end

% Calculate parameters
output.a=a; % Dirichlet parameters
output.r=a/sum(a); % (Expected) Model frequency
if Rec==1, output.RecIter=RecIter; end
% Exceedance Probability
EP=zeros(Nmodel,1);
for i=1:Nmodel
    EP(i)=integral(@(x)EPpdf(x,a(i),a(1:Nmodel~=i)),0,a(i))+...
        integral(@(x)EPpdf(x,a(i),a(1:Nmodel~=i)),a(i),Inf); % Integral of the posterior
end
output.EP=EP; % EP output

end

% Product of the gamcdf(s) for other models & the gampdf for the selected model
% Consistent with the algorithm in the Model Assessment, Comparison &
% Selection (MACS) Toolbox of SPM
% ## Reference ##
% - Soch, J. (2018), "Exceedance Probabilities for the Dirichlet Distribution", arXiv: 1611.01439
%
function p=EPpdf(x,a1,a2)
p=gampdf(x,a1,1); % Prob. of choosing a certain alpha
for i=1:length(a2)
    % Prob. that other alphas all <= the chosen alpha
    p=p.*gamcdf(x,a2(i),1);
end
end
