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
        if nargin==1
            Method.Algorithm='BinarySearch'; % binary search as default
            Method.Stat='CircSD';
            Method.Step=1;
        elseif isfield(Method,'Algorithm') && strcmp(Method.Algorithm,'BinarySearch')
            if ~isfield(Method,'Step')
                Method.Step=1;
            end
            if ~isfield(Method,'Stat')
                Method.Stat='circSD';
            end
        elseif isfield(Method,'Algorithm') && strcmp(Method.Algorithm,'KernelDensity')
            if ~isfield(Method,'Step')
                Method.Step=ceil((abs(Data.sample_range(end)-Data.sample_range(1)))/20);
            end
            if ~isfield(Method,'Stat')
                Method.Stat='CircSD';
            end
        end
        % Get curve
        f_raw=zeros(length(sample_range),1);
        for i=1:length(Data.sample_range)
            trial_id=Data.sample==Data.sample_range(i);
            errors=Data.response(trial_id)-Data.sample(trial_id);
            f_raw(i)=CircSummary_BMW(errors(trial_id),Method.Stat);
        end
        switch Method.Algorithm
            case 'BinarySearch'
                % Kernel average smoother
                kappa=100; % kernel radius
                f_smooth=zeros(length(Data.sample_range),1);
                for i=1:length(Data.sample_range)
                    f_weight=exp(kappa.*cosd(Data.sample_range-Data.sample_range(i)))./(2*pi*besseli(0,kappa));
                    f_weight=f_weight/sum(f_weight);
                    f_smooth(i)=sum(f_raw.*f_weight);
                end
            case 'KernelDensity'
            case 'FourierFit'
        end
end

end