%% Simulation_BMW
%
% Generate data based on a given likelihood function (e.g. a cognitive model)
% -----------------------
% Data1=Simulation_BMW(Config,Model,Param,Data0,Data_List,Data_Range,Input)
%
% ## Input ##
%

function Data1=Simulation_BMW(Config,Model,Param,Data0,Data_List,Data_Range,Input,Fix_List)

%% Prologue
if nargin<6 || nargin>8
    error('The number of input variables is not valid...')
elseif ~isempty(Input) && ~isfield(Input,'Output')
    Input.Output='LLH';
end
if ~iscell(Data_List) || ~iscell(Data_Range)
    error('The name of the branch of the data and the range of each branch should be cell variables...')
else
    uniq_Data_List=unique(Data_List);
    % uniq_Data_Range=unique(Data_Range);
    Data_Length=zeros(length(uniq_Data_List),1);
    Data_Ind=[];
    for i=1:length(uniq_Data_List)
        Data_Length(i)=sum(strcmp(Data_List,uniq_Data_List{i}));
        Data_Ind=strcat(Data_Ind,i*ones(1,Data_Length(i)));
    end
end
% fixed data fields
fix=0;
if nargin==8
    fix=1;
    if ~iscell(Fix_List)
        error('sorry, the input is not valid...')
    else
        Nfix=length(Fix_List);
        for i=1:Nfix
            DataFix.(Fix_List{i})=Data0.(Fix_List{i});
        end
    end
end
if ~ischar(Model)
    error('The name of the model should be a string variable...')
end
if isfield(Config,'Transform')
    Transform=config.Transform;
else
    Transform='Probit';
end
% start values
if isfield(Config,'start')
    start=zeros(length(Data_List),1);
    discrete_list=zeros(length(Data_List),1);
    for i=1:length(Data_List)
        if iscell(Data_Range{i})
            start(i)=MCMCConvert_BMW(Config.start(i),max(Data_Range{i}),min(Data_Range{i}),Transform);
        else
            discrete_list(i)=1;
            ind_start=find(Data_Range{i}==Config.start(i));
            start(i)=MCMCConvert_BMW(ind_start,length(Data_Range{i}),0,Transform);
        end
    end
else
    start=zeros(length(Data_List),1);
    discrete_list=zeros(length(Data_List),1);
    for i=1:length(Data_List)
        if iscell(Data_Range{i})
            start(i)=0;
        else
            discrete_list(i)=1;
            ind_temp=0:length(Data_Range{i});
            start(i)=MCMCConvert_BMW(randsample(ind_temp,1),max(ind_temp),0,Transform);
        end
    end
end
if sum(discrete_list)==0
    discrete=0;
else
    discrete=1;
end
if isstruct(Data0)
    Data1=Data0;
    if fix==1
        Data0_Field=setdiff(fieldnames(Data0),Fix_List);
    else
        Data0_Field=fieldnames(Data0);
    end
end
Ntrial=Config.Ntrial;

%% Simulation
% We use MH-MCMC for data generation
for i=1:length(uniq_Data_List)
    Data1.(uniq_Data_List{i})=zeros(Ntrial,Data_Length(i));
end
NewFunction=str2func(Model);
state0_0=start;
% De-interpolation
if discrete==1
    state0=zeros(size(state0_0,1),size(state0_0,2));
    for i=1:length(Data_List)
        if discrete_list(i)==1
            state0_2=MCMCConvert_BMW(state0_0(i),length(Data_Range{i}),0,['Inverse',Transform]);
            [~, state_hahaha]=min(state0_2-1:length(Data_Range{i}));
            state0_ind=state_hahaha(1);
            state0(i)=Data_Range{i}(state0_ind);
        else
            state0(i)=state0_0(i);
        end
    end
else
    state0=state0_0;
end
% pre-determined fields
if fix==1
    data0=DataFix;
    data1=DataFix;
end
for i=1:length(Data0_Field)
    data0.(Data0_Field{i})=Data0.(Data0_Field{i})(1,:);
end
% fields to simulate
for i=1:length(uniq_Data_List)
    data0.(uniq_Data_List{i})=[];
end
% Set tuning parameters
if isfield(Config,'cov')
    cov=Config.cov;
else
    cov=eye(length(Data_List))*(2.38^2/length(Data_List));
end
% Start chain
for state=1:Ntrial
        % Transition
        state1_0=mvnrnd(state0_0,cov);
        % De-interpolation
        if discrete==1
            state1=zeros(size(state1_0,1),size(state1_0,2));
            for i=1:length(Data_List)
                if discrete_list(i)==1
                    state1_2=MCMCConvert_BMW(state1_0(i),length(Data_Range{i}),0,['Inverse',Transform]);
                    [~, state_hahaha]=min(state1_2-1:length(Data_Range{i}));
                    state1_ind=state_hahaha(1);
                    state1(i)=Data_Range{i}(state1_ind);
                else
                    state1(i)=state1_0(i);
                end
            end
        else
            state1=state1_0;
        end
        % Write data1 & data0
        % pre-determined fields
        for i=1:length(Data0_Field)
            data1.(Data0_Field{i})=Data0.(Data0_Field{i})(1:state,:);
            data0.(Data0_Field{i})=Data0.(Data0_Field{i})(1:state,:);
        end
        % fields to simulate
        for i=1:length(uniq_Data_List)
            data1.(uniq_Data_List{i})=vertcat(data0.(uniq_Data_List{i}),state1(Data_Ind==i));
            data0.(uniq_Data_List{i})=vertcat(data0.(uniq_Data_List{i}),state0(Data_Ind==i));
        end
        % Metropolis criterion
        minMHr=rand; % sample Metropolis-Hastings ratio threshold
        % Get likelihood
        PropTar=NewFunction(Param, data1, Input);
        PrevTar=NewFunction(Param, data0, Input);
        MHr=exp(-PropTar+PrevTar)*(mvnpdf(state0_0,state1_0,cov)/mvnpdf(state1_0,state0_0,cov)); % current MH ratio
        if MHr>minMHr
            fprintf('accept!\n\n\n\nhahahaahahahha\n\n')
            state0_0=state1_0; % accept proposal
            state0=state1;
            data0=data1;
        end
        state0_0
        state0
        % Write Data1
        for i=1:length(uniq_Data_List)
            Data1.(uniq_Data_List{i})(state)=data0.(uniq_Data_List{i})(state);
        end
        state
end

end
