%% Categorization (For Delayed Recall Task)
%
%

function Data=Categorization_BMW(Data,Method)

if ~isstruct(Data)
    error('The input data is not valid...')
end
if isfield(Data,'category_range')
    if isfield(Data,'sample')
        LearnType='Supervised';
    else
        error('The input data is not valid...')
    end
else
    if isfield(Data,'sample') && isfield(Data,'response')
        LearnType='Unsupervised';
    else
        error('The input data is not valid...')
    end
end
nt=0; % non-target
if isfield(Data,'sample_nt')
    nt=1;
    Nnt=size(Data.sample_nt,2);
end

Ntrial=length(Data.sample);
Data.category=zeros(Ntrial,1);
switch LearnType
    case 'Supervised'
        cat_diff0=abs(CircDist_BMW('Diff',repmat(Data.category_range,[length(Data.sample),1]),repmat(Data.sample,[1,length(Data.category_range)])));
        [~,cat_ind0]=min(cat_diff0,[],2);
        Data.category=Data.category_range(cat_ind0)';
        if nt==1
            Data.category_nt=zeros(Ntrial,Nnt);
            for i=1:Nnt
                cat_diff0=abs(CircDist_BMW('Diff',repmat(Data.category_range,[length(Data.sample_nt(:,i)),1]),repmat(Data.sample_nt(:,i),[1,length(Data.category_range)])));
                [~,cat_ind0]=min(cat_diff0,[],2);
                Data.category_nt(:,i)=Data.category_range(cat_ind0)';
            end
        end
    case 'Unsupervised'
        
end

end